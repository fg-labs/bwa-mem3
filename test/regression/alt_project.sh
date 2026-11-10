#!/usr/bin/env bash
# test/regression/alt_project.sh
#
# End-to-end behavior of --alt-project.
#
# ALT-aware mapping keeps a primary-assembly hit primary whatever an ALT contig
# scores, so a read whose haplotype exists only as an ALT is placed by the
# primary assembly alone. --alt-project lifts ALT hits onto the primary assembly
# through <idxbase>.alt before pairing and MAPQ, so the ALT decides.
#
# Fixture (generated): copy X on chrA and a paralog Y on chrB differing at 8
# sites. chrA_alt1 is a haplotype of X that carries Y's base at 6 of them (gene
# conversion), one private SNP and a 3-bp deletion. Read pairs come from
# chrA_alt1. Reads covering the private SNP match the ALT best; on the primary
# assembly alone they cannot tell X from Y.
#
# Asserts:
#   * without the flag, two runs give identical output (determinism; flag-off
#     byte-identity to earlier releases is pinned by the rest of the suite), and
#     the SNP reads are not all at chrA with MAPQ > 0 (so the test is not vacuous)
#   * with --alt-project (paired, -I fixed): every SNP-covering read is on chrA
#     at its lifted position (chrA:3001 + its ALT start, from the read name) as
#     150M with MAPQ > 0, both mates of its pair are on chrA as a proper pair, its
#     AS is the ALT score that placed it (150), and its NM is exactly the number
#     of bases where its ALT window differs from chrA; no record (RNAME, RNEXT or
#     any tag) names the ALT contig; and every read spanning the ALT deletion
#     carries the 3D at the lifted position (the read name encodes its ALT start,
#     an oracle independent of the aligner)
#   * a reverse-complemented ALT (.alt FLAG 16) gives the same records
#   * ALTs whose insertion and deletion cancel inside a read (1I...1D, 2I...2D,
#     5I...5D, 40 bp apart; 35I...35D, 30 bp apart, the widest excursion the
#     re-alignment band of a 150 bp read can follow) and one with a wide pair
#     (25I then 30D): every lifted read is written at its lifted position with
#     both indels and NM equal to the indel bases, not as an ungapped run of
#     mismatches; for the 1-5 bp balanced pairs its record equals the default
#     run's chrA record apart from AS and the tags that name the ALT hit
#   * a 40I...40D ALT, past that band: its reads stay on the ALT, with records
#     equal to the default run's, not lifted with a garbled CIGAR
#   * a read that runs off the end of a divergent ALT onto sequence matching
#     chrA: its shorter ALT hit outscores its full chrA hit, but not by more than
#     the 50-54 chrA-aligned bases it leaves out are worth, so the full chrA record
#     is kept, unclipped: equal to the default run's primary record apart from
#     the SA and pa tags, with no ALT supplementary record (the default writes
#     one, so the check is not vacuous)
#   * the other direction: a read whose ALT hit misses 4 bases at its start
#     (4S146M, AS 146) but beats its full chrA hit (AS 120) by more than they are
#     worth is placed there by the lifted hit, not at a paralog elsewhere on chrA
#     (AS 130) where the default places it
#   * a tandem repeat U U U' on chrA whose ALT is U U' U': a read across copies 1
#     and 2 matches the ALT at its own locus and chrA one unit on equally well, so
#     both placements survive, as MAPQ 0 with XS equal to AS and, under -a, a
#     secondary record one unit away
#   * a fixture whose ALT starts with 300 bp of ALT-only sequence (.alt CIGAR
#     300S1000M): reads inside that clip stay on the ALT contig; reads straddling
#     its edge keep their ALT hit on the ALT (records equal to the default run's,
#     no long insertion); reads past it are at chrA or chrB at the lifted position
#     and no record of theirs names the ALT, which the default run's do (so the
#     check is not vacuous); an unmapped read fails the check
#   * records are byte-identical (apart from @PG) at -t 1, 2 and 4; --fast at
#     -t 1 and 4, where it also lifts the SNP reads to their chrA positions and
#     no record names the ALT; and --dedup-reads on and off over input with every
#     pair duplicated, with the memo shown to engage (dup > 0, final=ON) and the
#     BWAMEM3_DEDUP_READS_VERIFY check of duplicate regions passing
#   * --alt-project is refused with -j, with --meth, without a .alt and with a
#     .alt that has no usable line, and --compat-allow-divergent admits it under
#     --compat (the --compat refusal itself is checked in compat_byte_identical.sh)
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
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

python3 - << 'EOF'
import random
random.seed(7)
B = "ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(c): return random.choice([x for x in B if x != c])
RC = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(RC)[::-1]
X = rnd(2000)
chrA = rnd(3000) + X + rnd(3000)
psv = sorted(random.sample(range(200, 1800), 8))
Y = list(X)
for p in psv: Y[p] = mut(X[p])
chrB = rnd(3000) + "".join(Y) + rnd(3000)
alt = list(X)
for p in psv[:6]: alt[p] = Y[p]
alt[1000] = mut(X[1000])                        # private SNP
# ALT positions (all before the deletion, so equal to X positions) where the ALT
# differs from chrA: the oracle for a lifted read's NM.
open("diffs.txt", "w").write(" ".join(str(p) for p in sorted(psv[:6] + [1000])) + "\n")
alt = "".join(alt[:1500] + alt[1503:])          # 3-bp deletion at X:1500
cig = "1500M3D497M"
for name, seq, flag in (("fwd", alt, 0), ("rev", rc(alt), 16)):
    open(f"{name}.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>chrA_alt1\n{seq}\n")
    open(f"{name}.fa.alt", "w").write(f"chrA_alt1\t{flag}\tchrA\t3001\t60\t{cig}\t*\t0\t0\t*\t*\n")
open("noalt.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>chrA_alt1\n{alt}\n")
with open("r1.fq", "w") as r1, open("r2.fq", "w") as r2:
    for i in range(60):
        s = random.randint(0, len(alt) - 400)
        f = alt[s:s + 400]
        q = "I" * 150
        r1.write(f"@p{i}_{s}\n{f[:150]}\n+\n{q}\n")
        r2.write(f"@p{i}_{s}\n{rc(f[-150:])}\n+\n{q}\n")
# Clip fixture: 300 bp of ALT-only sequence, then X[0:1000] at chrA:3001.
alt2 = rnd(300) + X[:1000]
open("clip.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>chrA_alt2\n{alt2}\n")
open("clip.fa.alt", "w").write("chrA_alt2\t0\tchrA\t3001\t60\t300S1000M\t*\t0\t0\t*\t*\n")
with open("clip.fq", "w") as fq:
    for s in list(range(0, 300, 10)) + list(range(400, 1100, 50)):
        fq.write(f"@c_{s}\n{alt2[s:s + 150]}\n+\n{'I' * 150}\n")
# Indel fixtures: an ALT of Z = rnd(2000) at chrA:3001 with an insertion after
# Z[1000) and a deletion GAP bases later. Single-end reads from the ALT start at
# s < 1000 and cover both indels; their lifted position is chrA:3001 + s.
def indel_fixture(name, ins, gap, dl, seed):
    rng = random.Random(seed)
    r = lambda n: "".join(rng.choice(B) for _ in range(n))
    Z = r(2000)
    a = Z[:1000] + r(ins) + Z[1000:1000 + gap] + Z[1000 + gap + dl:]
    cig = f"1000M{ins}I{gap}M{dl}D{2000 - 1000 - gap - dl}M"
    open(f"{name}.fa", "w").write(f">chrA\n{r(3000) + Z + r(3000)}\n>chrA_alt1\n{a}\n")
    open(f"{name}.fa.alt", "w").write(f"chrA_alt1\t0\tchrA\t3001\t60\t{cig}\t*\t0\t0\t*\t*\n")
    with open(f"{name}.fq", "w") as fq:
        for s in (930, 950, 970):
            fq.write(f"@i_{s}_{ins}_{dl}\n{a[s:s + 150]}\n+\n{'I' * 150}\n")
for k in (1, 2, 5):
    indel_fixture(f"bal{k}", k, 40, k, 100 + k)
indel_fixture("bal35", 35, 30, 35, 135)
indel_fixture("bal40", 40, 30, 40, 140)
indel_fixture("wide", 25, 30, 30, 12)
# End fixture: an ALT of W[0:1000) at chrA:3001 with 12 ALT-only SNPs in its last
# 100 bases. Each read is the ALT's last 100 - k bases then W's next 50 + k, which
# match chrA exactly: 100 - k on the ALT (no mismatch) against 150 on chrA (12
# mismatches).
# Its own generator, so the fixtures after it are as before.
rng = random.Random(11)
rw = lambda n: "".join(rng.choice(B) for _ in range(n))
W = rw(2000)
altw = list(W[:1000])
for p in range(904, 1000, 8): altw[p] = rng.choice([x for x in B if x != altw[p]])
altw = "".join(altw)
open("end.fa", "w").write(f">chrA\n{rw(3000) + W + rw(3000)}\n>chrA_alt1\n{altw}\n")
open("end.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1000M\t*\t0\t0\t*\t*\n")
with open("end.fq", "w") as fq:
    for k in (0, 2, 4):
        seq = altw[900 + k:1000] + W[1000:1050 + k]
        fq.write(f"@e_{900 + k}\n{seq}\n+\n{'I' * len(seq)}\n")
# Duplicated pairs (new names) for the --dedup-reads check.
with open("r1d.fq", "w") as o1, open("r2d.fq", "w") as o2:
    for src, dst in (("r1.fq", o1), ("r2.fq", o2)):
        lines = open(src).read().splitlines()
        for i in range(0, len(lines), 4):
            dst.write("\n".join(lines[i:i + 4]) + "\n")
            dst.write("\n".join([lines[i] + "_dup"] + lines[i + 1:i + 4]) + "\n")
# Tandem fixture: U U U' on chrA, U U' U' on the ALT (U' = U with one SNP).
U = rnd(100)
U2 = U[:20] + mut(U[20]) + U[21:]
lf, rt = rnd(500), rnd(500)
open("tan.fa", "w").write(f">chrA\n{rnd(3000) + lf + U + U + U2 + rt + rnd(3000)}\n>chrA_alt1\n{lf + U + U2 + U2 + rt}\n")
open("tan.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1300M\t*\t0\t0\t*\t*\n")
alt3 = lf + U + U2 + U2 + rt
with open("tan.fq", "w") as fq:
    for s in (500, 510, 530, 550):  # across copies 1 and 2, over U2's SNP
        fq.write(f"@t_{s}\n{alt3[s:s + 150]}\n+\n{'I' * 150}\n")
# Paralog fixture (its own generator): an ALT of V[0:1000) at chrA:3001 with
# SNPs at s+1, s+3 and six in s+104..s+144, s = 400. The read is V[s:s+100] then
# the ALT's next 50 bases: on the ALT it aligns from base 4 (4S146M, AS 146), on
# chrA's copy of V over all 150 bases with the six SNPs (AS 120), and on a
# paralog P (the read with 4 mismatches, at chrA:8001) at AS 130.
rng = random.Random(5)
rv = lambda n: "".join(rng.choice(B) for _ in range(n))
mutv = lambda c: rng.choice([x for x in B if x != c])
V = rv(2000)
s = 400
altv = list(V[:1000])
for p in [s + 1, s + 3] + list(range(s + 104, s + 150, 8)): altv[p] = mutv(altv[p])
altv = "".join(altv)
read = V[s:s + 100] + altv[s + 100:s + 150]
P = list(read)
for p in (20, 60, 90, 130): P[p] = mutv(P[p])
open("par.fa", "w").write(f">chrA\n{rv(3000) + V + rv(3000) + ''.join(P) + rv(1000)}\n>chrA_alt1\n{altv}\n")
open("par.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1000M\t*\t0\t0\t*\t*\n")
open("par.fq", "w").write(f"@v_{s}\n{read}\n+\n{'I' * 150}\n")
# The reader every check below shares: a SAM's records in file order, grouped by
# QNAME, the primary record of each read, and a record's tags.
open("samrecs.py", "w").write('''def records(path):
    for l in open(path):
        if l[0] != "@": yield l.rstrip("\\n").split("\\t")
def recs(path):
    out = {}
    for t in records(path): out.setdefault(t[0], []).append(t)
    return out
def prim(path):
    return {t[0]: t for t in records(path) if not int(t[1]) & 0x900}
def tags(t):
    return {x[:2]: x[5:] for x in t[11:]}
''')
EOF

for ref in fwd rev noalt clip bal1 bal2 bal5 bal35 bal40 wide tan end par; do
    "$BIN" index "$ref.fa" > "idx.$ref.log" 2>&1 || fail "index $ref failed"
done
# An index whose .alt has no usable line (its CIGAR length does not match the ALT).
for f in fwd.fa.*; do
    case "$f" in
        *.alt) ;;
        *) ln -s "$f" "bad.fa${f#fwd.fa}" ;;
    esac
done
ln -s fwd.fa bad.fa
printf 'chrA_alt1\t0\tchrA\t3001\t60\t100M\t*\t0\t0\t*\t*\n' > bad.fa.alt

aln() { # ref out extra...
    local ref=$1 out=$2
    shift 2
    "$BIN" mem -t 2 -I 400,10 "$@" "$ref.fa" r1.fq r2.fq > "$out" 2> "$out.log" || fail "mem $* on $ref failed (see $out.log)"
}
# Summarize primary records: name, mate, chrom, pos, MAPQ, CIGAR, FLAG, AS, NM,
# and the mate's ALT span (R1 = [s, s+150), R2 = [s+250, s+400)).
summ() {
    python3 - "$1" << 'EOF'
import sys
from samrecs import records, tags
for t in records(sys.argv[1]):
    flag = int(t[1])
    if flag & 0x900: continue
    s = int(t[0].split("_")[1]); m = 1 if flag & 0x40 else 2
    a, b = (s, s + 150) if m == 1 else (s + 250, s + 400)
    tg = tags(t)
    print(t[0], m, t[2], t[3], t[4], t[5], flag, tg.get("AS", "NA"), tg.get("NM", "NA"), a, b)
EOF
}

aln fwd def1.sam
aln fwd def2.sam
grep -v '^@PG' def1.sam > def1.nopg
grep -v '^@PG' def2.sam > def2.nopg
cmp -s def1.nopg def2.nopg || fail "default output is not deterministic"
grep -q "alt-project" def1.sam.log && fail "default run reports alt-project"

aln fwd proj.sam --alt-project
grep -q "alt-project: 1 ALT contigs" proj.sam.log || fail "--alt-project did not report the lift table"
aln rev projr.sam --alt-project

summ def1.sam > def.tsv
summ proj.sam > proj.tsv
summ projr.sam > projr.tsv

# SNP-covering mates (ALT coordinate 1000).
n_snp=$(awk '$10 <= 1000 && 1000 < $11' proj.tsv | wc -l)
[[ $n_snp -ge 5 ]] || fail "fixture has only $((n_snp)) SNP-covering reads"
bad=$(awk '$10 <= 1000 && 1000 < $11 && ($3 != "chrA" || $5 == 0)' proj.tsv | wc -l)
[[ $bad -eq 0 ]] || fail "--alt-project left $((bad)) SNP-covering reads off chrA or at MAPQ 0"
def_ok=$(awk '$10 <= 1000 && 1000 < $11 && $3 == "chrA" && $5 > 0' def.tsv | wc -l)
[[ $def_ok -lt $n_snp ]] || fail "default already places every SNP read at chrA with MAPQ > 0 (test is vacuous)"
# AS is the ALT score that placed the read; the read is 150M at its lifted
# position (chrA:3001 + its ALT start), and NM is exactly the number of bases
# where its ALT window differs from chrA.
cat > snpcheck.py << 'EOF'
import sys
diffs = [int(x) for x in open(sys.argv[2]).read().split()]
n = 0
for l in open(sys.argv[1]):
    name, m, c, pos, q, cig, flag, AS, NM, a, b = l.split()
    a, b = int(a), int(b)
    if not (a <= 1000 < b): continue
    n += 1
    want_nm = sum(a <= d < b for d in diffs)
    if c != "chrA" or int(pos) != 3001 + a or cig != "150M" or AS != "150" or int(NM) != want_nm:
        print(f"{name}/{m}: {c}:{pos} {cig} AS:{AS} NM:{NM}, want chrA:{3001 + a} 150M AS:150 NM:{want_nm}",
              file=sys.stderr)
        sys.exit(1)
if n < 5: print(f"only {n} SNP-covering reads", file=sys.stderr); sys.exit(1)
EOF
python3 snpcheck.py proj.tsv diffs.txt || fail "a lifted SNP read is misplaced or has the wrong AS or NM"
# Both mates of every SNP fragment: on chrA, proper pair.
awk '$10 <= 1000 && 1000 < $11 { print $1 }' proj.tsv | sort -u > snp.names
bad=$(awk 'NR == FNR { n[$1] = 1; next } ($1 in n) && ($3 != "chrA" || int($7 / 2) % 2 != 1)' snp.names proj.tsv | wc -l)
[[ $bad -eq 0 ]] || fail "$((bad)) mates of SNP fragments are off chrA or not a proper pair"
# No record names the ALT contig, in any field.
grep -v '^@' proj.sam > proj.body
grep -q 'chrA_alt1' proj.body && fail "--alt-project wrote a record that names the ALT contig"
# Deletion: X:1500 is chrA 0-based 4500, so the D starts at 1-based 4501. Every
# mate with at least 10 ALT bases on each side of it must carry it there.
python3 - proj.tsv << 'EOF' || fail "a read spanning the ALT deletion lacks the 3D at chrA:4501"
import re, sys
n = 0
for l in open(sys.argv[1]):
    name, m, c, pos, q, cig, flag, AS, NM, a, b = l.split()
    if not (int(a) <= 1490 and int(b) >= 1510): continue
    n += 1
    ref, ok = 0, False
    for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        if op == "D" and int(ln) == 3 and int(pos) + ref == 4501: ok = True
        if op in "MDN=X": ref += int(ln)
    if c != "chrA" or not ok:
        print(f"{name}/{m}: {c}:{pos} {cig}", file=sys.stderr); sys.exit(1)
if n < 3: print(f"only {n} reads span the deletion", file=sys.stderr); sys.exit(1)
EOF
cmp -s proj.tsv projr.tsv || fail "reverse-complemented ALT (.alt FLAG 16) placed reads differently"

# Indel fixtures (single-end): every lifted read keeps both ALT indels at its
# lifted position, with NM equal to the indel bases (the reads are exact ALT
# substrings), rather than an ungapped run of mismatches. For the 1-5 bp
# balanced pairs the record also equals the default run's chrA record, apart
# from AS (the ALT score) and the tags that name the ALT hit (at 35 bp the
# default run soft-clips the read on chrA instead of paying for two gaps).
for ref in bal1 bal2 bal5 bal35 wide; do
    "$BIN" mem -t 2 --alt-project "$ref.fa" "$ref.fq" > "$ref.proj.sam" 2> "$ref.proj.log" \
        || fail "mem --alt-project on $ref failed"
    "$BIN" mem -t 2 "$ref.fa" "$ref.fq" > "$ref.def.sam" 2> "$ref.def.log" || fail "mem on $ref failed"
    python3 - "$ref.proj.sam" "$ref.def.sam" "$ref" << 'EOF' || fail "--alt-project misaligned a read across the ALT indels of $ref"
import re, sys
from samrecs import prim, tags
proj, dflt = prim(sys.argv[1]), prim(sys.argv[2])
if len(proj) != 3: print(f"{len(proj)} primary records, want 3", file=sys.stderr); sys.exit(1)
for name, t in proj.items():
    s, ins, dl = (int(x) for x in name.split("_")[1:])
    tg = tags(t)
    ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
    ok = (t[2] == "chrA" and int(t[3]) == 3001 + s and (str(ins), "I") in ops and (str(dl), "D") in ops
          and all(op in "MID" for _, op in ops) and int(tg.get("NM", -1)) == ins + dl and tg.get("AS") == "150")
    if ok and ins == dl and ins <= 5:
        d = dflt[name]
        dtg = tags(d)
        ok = d[2:6] == t[2:6] and dtg.get("MD") == tg.get("MD") and dtg.get("NM") == tg.get("NM")
    if not ok:
        print(f"{name}: {t[2]}:{t[3]} {t[5]} NM:{tg.get('NM')} AS:{tg.get('AS')}", file=sys.stderr)
        sys.exit(1)
EOF
done

# Past the band: a 40I...40D ALT needs a band of 40, more than the 35 the
# re-alignment of a 150 bp read can have, so its reads stay on the ALT and their
# records equal the default run's.
"$BIN" mem -t 2 --alt-project bal40.fa bal40.fq > bal40.proj.sam 2> bal40.proj.log || fail "mem --alt-project on bal40 failed"
"$BIN" mem -t 2 bal40.fa bal40.fq > bal40.def.sam 2> bal40.def.log || fail "mem on bal40 failed"
grep -v '^@PG' bal40.proj.sam > bal40.proj.nopg
grep -v '^@PG' bal40.def.sam > bal40.def.nopg
cmp -s bal40.proj.nopg bal40.def.nopg || fail "--alt-project changed reads whose ALT indels are past the re-alignment band (bal40)"
grep -q $'\tchrA_alt1\t' bal40.def.nopg || fail "bal40: the default run has no ALT record (test is vacuous)"

# End fixture (single-end): for read e_<900 + k> (k = 0, 2, 4) the ALT hit is
# 100 - k bases at AS 100 - k and the chrA hit 150M at AS 90 + k (its first
# mismatch is at read offset 4 - k; AS is the local score, which leaves out the
# bases up to it, while the clipping penalty keeps the CIGAR end to end); the
# 50 + k bases the lifted hit leaves out are worth 50 + k, so it does not
# replace the chrA hit and is dropped. The read's only record is the default
# run's primary chrA record (150M at chrA:3001 + its ALT start),
# field for field and tag for tag, minus what the default derives from the ALT
# hit: its SA tag (to the default's ALT supplementary record, which is gone too)
# and its pa tag. The default run must have that supplementary record and those
# tags, so the comparison is not vacuous.
"$BIN" mem -t 2 --alt-project end.fa end.fq > end.proj.sam 2> end.proj.log || fail "mem --alt-project on end failed"
"$BIN" mem -t 2 end.fa end.fq > end.def.sam 2> end.def.log || fail "mem on end failed"
python3 - end.proj.sam end.def.sam << 'EOF' || fail "--alt-project let a shorter lifted hit replace a full chrA hit"
import sys
from samrecs import recs, tags
got, dflt = recs(sys.argv[1]), recs(sys.argv[2])
if len(got) != 3 or set(got) != set(dflt): sys.exit(f"{len(got)} reads, want the default's 3")
for name, rs in got.items():
    s = int(name.split("_")[1])
    dp = [t for t in dflt[name] if not int(t[1]) & 0x900]
    ds = [t for t in dflt[name] if int(t[1]) & 0x800 and t[2] == "chrA_alt1"]
    if len(dp) != 1 or len(ds) != 1: sys.exit(f"{name}: the default run has no ALT supplementary record (vacuous)")
    d = dp[0]
    if "SA" not in tags(d) or "pa" not in tags(d):
        sys.exit(f"{name}: the default primary record has no SA or pa tag (vacuous)")
    want = d[:11] + [x for x in d[11:] if not x.startswith(("SA:Z:", "pa:f:"))]
    if len(rs) != 1 or rs[0] != want or rs[0][2] != "chrA" or int(rs[0][3]) != 3001 + s or rs[0][5] != "150M":
        got_s = ["\t".join(t[1:6] + t[11:]) for t in rs]
        sys.exit(f"{name}: {got_s}, want one record, the default primary minus SA and pa: {want[1:6] + want[11:]}")
EOF

# Paralog fixture (single-end), the other direction: the lifted hit misses the
# read's first 4 bases (4S146M, AS 146) and the chrA hit at the same locus aligns
# all 150 (AS 120); 146 beats 120 + 4, so the lifted hit replaces it, and the read
# is written at chrA:3405 (3001 + 400 + 4) rather than at the paralog (chrA:8001,
# AS 130), where the default run places it.
"$BIN" mem -t 2 --alt-project par.fa par.fq > par.proj.sam 2> par.proj.log || fail "mem --alt-project on par failed"
"$BIN" mem -t 2 par.fa par.fq > par.def.sam 2> par.def.log || fail "mem on par failed"
python3 - par.proj.sam par.def.sam << 'EOF' || fail "--alt-project did not let a much better lifted hit that clips 4 bases replace a full chrA hit"
import sys
from samrecs import recs, prim, tags
got, d = recs(sys.argv[1]), prim(sys.argv[2])
if set(got) != {"v_400"}: sys.exit(f"reads {sorted(got)}, want v_400")
rs = got["v_400"]
p = [t for t in rs if not int(t[1]) & 0x900]
if not p: sys.exit("no primary record")
p = p[0]
if (p[2], p[3], p[5], tags(p).get("AS")) != ("chrA", "3405", "4S146M", "146") or int(p[4]) == 0:
    sys.exit(f"primary {p[2]}:{p[3]} {p[5]} MAPQ {p[4]} AS:{tags(p).get('AS')}, want chrA:3405 4S146M AS:146, MAPQ > 0")
if any("chrA_alt1" in "\t".join(t) for t in rs): sys.exit("a record names the ALT")
dv = d["v_400"]
if (dv[2], dv[3]) != ("chrA", "8001"): sys.exit(f"the default places it at {dv[2]}:{dv[3]}, not the paralog (vacuous)")
EOF

# Tandem fixture (single-end): the lifted hit and the chrA placement one repeat
# unit away score the same, so both survive: MAPQ 0, XS equal to AS, and under
# -a a secondary at the other placement.
"$BIN" mem -t 2 -a --alt-project tan.fa tan.fq > tan.sam 2> tan.log || fail "mem -a --alt-project on tan failed"
python3 - tan.sam << 'EOF' || fail "--alt-project dropped a competing placement in a tandem repeat"
import sys
from samrecs import recs, tags
by_name = recs(sys.argv[1])
for rs in by_name.values():
    for t in rs:
        if int(t[1]) & 0x4: print(f"{t[0]}: unmapped", file=sys.stderr); sys.exit(1)
if len(by_name) != 4: print(f"{len(by_name)} reads, want 4", file=sys.stderr); sys.exit(1)
for name, rs in by_name.items():
    s = int(name.split("_")[1])
    p = [t for t in rs if not int(t[1]) & 0x900]
    if not p: print(f"{name}: no primary record", file=sys.stderr); sys.exit(1)
    tg = tags(p[0])
    pos150 = {(t[2], int(t[3])) for t in rs if not int(t[1]) & 0x800 and "AS:i:150" in t[11:]}
    if (len(p) != 1 or p[0][4] != "0" or tg.get("AS") != "150" or tg.get("XS") != "150"
            or pos150 != {("chrA", 3001 + s), ("chrA", 3101 + s)}):
        print(f"{name}: MAPQ {p[0][4]} AS:{tg.get('AS')} XS:{tg.get('XS')} at {sorted(pos150)}", file=sys.stderr)
        sys.exit(1)
EOF

# Clip fixture (single-end). Reads wholly inside the ALT-only 300 bp stay on the
# ALT. Reads straddling its edge keep their ALT hit on the ALT: their records
# equal the default run's and none carries a long insertion (an unguarded lift
# writes the ALT-only bases as one). Reads past it are at chrA or chrB at the
# lifted position, and none of their records names the ALT: the ALT hit went to
# chrA and merged with chrA's own. The default run names the ALT on those reads
# (as XA), so that last check fails there and is not vacuous. An unmapped read is
# a failure.
"$BIN" mem -t 2 --alt-project clip.fa clip.fq > clip.sam 2> clip.log || fail "mem --alt-project on clip failed"
"$BIN" mem -t 2 clip.fa clip.fq > clip.def.sam 2> clip.def.log || fail "mem on clip failed"
cat > clipcheck.py << 'EOF'
import re, sys
from samrecs import recs
got, dflt = recs(sys.argv[1]), recs(sys.argv[2])
n_in = n_edge = n_past = 0
for name, rs in got.items():
    s = int(name.split("_")[1])
    for t in rs:
        if int(t[1]) & 0x4: sys.exit(f"{name}: unmapped")
        if max([int(x) for x in re.findall(r"(\d+)I", t[5])] + [0]) >= 20: sys.exit(f"{name}: {t[5]}")
    p = [t for t in rs if not int(t[1]) & 0x900]
    if not p: sys.exit(f"{name}: no primary record")
    p = p[0]
    if s + 150 <= 300:
        n_in += 1
        if p[2] != "chrA_alt2": sys.exit(f"{name}: {p[2]}, want chrA_alt2")
    elif s < 300:
        n_edge += 1
        if not any(t[2] == "chrA_alt2" for t in rs): sys.exit(f"{name}: its ALT hit left the ALT")
        if rs != dflt[name]: sys.exit(f"{name}: records differ from the default run's")
    else:
        n_past += 1
        if p[2] not in ("chrA", "chrB") or int(p[3]) != 3001 + s - 300: sys.exit(f"{name}: {p[2]}:{p[3]}")
        if any("chrA_alt2" in "\t".join(t) for t in rs): sys.exit(f"{name}: a record names chrA_alt2")
if n_in < 10 or n_edge < 10 or n_past < 10: sys.exit(f"only {n_in} clip / {n_edge} edge / {n_past} past reads")
EOF
python3 clipcheck.py clip.sam clip.def.sam || fail "--alt-project mishandled reads in, across or past an ALT clip"
# The same check on the default run must fail, and for the reason it exists:
# a past-the-clip read whose records name the ALT.
python3 clipcheck.py clip.def.sam clip.def.sam 2> clipvac.log && fail "the clip check passes on the default run (vacuous)"
grep -q "a record names chrA_alt2" clipvac.log || fail "the clip check fails on the default run for another reason (see clipvac.log)"

# Determinism: records are byte-identical (apart from @PG, whose CL records the
# arguments) across thread counts, under --fast (its pdqsort orders score ties
# differently), and with --dedup-reads on and off over input with every pair
# duplicated, so memoization replays records.
nopg() { grep -v '^@PG' "$1"; }
for t in 1 4; do
    "$BIN" mem -t "$t" -I 400,10 --alt-project fwd.fa r1.fq r2.fq > "det.t$t.sam" 2> "det.t$t.log" || fail "mem -t $t --alt-project failed"
    "$BIN" mem -t "$t" -I 400,10 --fast --alt-project fwd.fa r1.fq r2.fq > "detf.t$t.sam" 2> "detf.t$t.log" || fail "mem -t $t --fast --alt-project failed"
    BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t "$t" -I 400,10 --dedup-reads on --alt-project fwd.fa r1d.fq r2d.fq \
        > "detd.on.t$t.sam" 2> "detd.on.t$t.log" || fail "mem -t $t --dedup-reads on --alt-project failed"
    "$BIN" mem -t "$t" -I 400,10 --dedup-reads off --alt-project fwd.fa r1d.fq r2d.fq > "detd.off.t$t.sam" 2> "detd.off.t$t.log" \
        || fail "mem -t $t --dedup-reads off --alt-project failed"
done
cmp -s <(nopg det.t1.sam) <(nopg det.t4.sam) || fail "--alt-project output differs between -t 1 and -t 4"
cmp -s <(nopg det.t1.sam) <(nopg proj.sam) || fail "--alt-project output differs between -t 1 and -t 2"
cmp -s <(nopg detf.t1.sam) <(nopg detf.t4.sam) || fail "--fast --alt-project output differs between -t 1 and -t 4"
# --fast lifts too: its SNP reads pass the same placement check as the default
# mode's (which also rules out an empty or unmapped output), and no record names
# the ALT (a file, not a pipe into grep -q, whose early exit could SIGPIPE the
# producer and, under pipefail, hide a match).
summ detf.t1.sam > detf.tsv
python3 snpcheck.py detf.tsv diffs.txt || fail "--fast --alt-project did not lift the SNP reads to their chrA positions"
grep -v '^@' detf.t1.sam > detf.body
grep -q 'chrA_alt1' detf.body && fail "--fast --alt-project wrote a record that names the ALT contig"
# --dedup-reads: byte-identical on and off (its documented contract), and the
# memo must engage, or on == off is vacuous: its stats line reports duplicate
# pairs and the final state ON.
for t in 1 4; do
    cmp -s <(nopg "detd.on.t$t.sam") <(nopg "detd.off.t$t.sam") || fail "--alt-project output differs with --dedup-reads on and off (-t $t)"
    grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* .*final=ON' "detd.on.t$t.log" \
        || fail "--dedup-reads on --alt-project did not memoize any pair (-t $t, see detd.on.t$t.log)"
done
cmp -s <(nopg detd.on.t1.sam) <(nopg detd.on.t4.sam) || fail "--dedup-reads on --alt-project output differs between -t 1 and -t 4"
# The VERIFY diagnostic aligns each duplicate itself and aborts if its regions
# (lifted ones and their band floor included) differ from its representative's;
# its output must equal the off run's, and it must have compared duplicates.
BWAMEM3_DEDUP_READS_VERIFY=1 BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t 4 -I 400,10 --dedup-reads on --alt-project fwd.fa r1d.fq r2d.fq \
    > detd.verify.sam 2> detd.verify.log || fail "--dedup-reads VERIFY found a duplicate whose regions differ (see detd.verify.log)"
grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* ' detd.verify.log \
    || fail "--dedup-reads VERIFY compared no duplicate pair (see detd.verify.log)"
cmp -s <(nopg detd.verify.sam) <(nopg detd.off.t4.sam) || fail "--dedup-reads VERIFY output differs from the off run's"

# Refusals
"$BIN" mem -j --alt-project fwd.fa r1.fq r2.fq > /dev/null 2> j.log && fail "-j --alt-project was accepted"
grep -q -e "cannot be combined with -j (it needs" j.log || fail "-j --alt-project gave the wrong error"
"$BIN" mem --meth --alt-project fwd.fa r1.fq r2.fq > /dev/null 2> m.log && fail "--meth --alt-project was accepted"
grep -q -e "cannot be combined with --meth" m.log || fail "--meth --alt-project gave the wrong error"
"$BIN" mem --alt-project noalt.fa r1.fq r2.fq > /dev/null 2> n.log && fail "--alt-project without a .alt was accepted"
grep -q "needs noalt.fa.alt" n.log || fail "--alt-project without a .alt gave the wrong error"
"$BIN" mem --alt-project bad.fa r1.fq r2.fq > /dev/null 2> b.log && fail "--alt-project with no usable .alt line was accepted"
grep -q "no ALT contig in bad.fa.alt has a usable alignment" b.log || fail "--alt-project with no usable .alt line gave the wrong error"
# The --compat refusal (both targets) is checked in compat_byte_identical.sh,
# with the other levers --compat refuses; only the override is checked here.
"$BIN" mem --compat bwa-mem --compat-allow-divergent --alt-project fwd.fa r1.fq r2.fq > /dev/null 2> ca.log \
    || fail "--compat-allow-divergent did not admit --alt-project (see ca.log)"

echo "PASS: alt_project ($((n_snp)) SNP reads on chrA as proper pairs, FLAG 16 identical, ALT indels kept across a lift, past-band indels kept on the ALT, full primary hits kept, a clearly better clipped lifted hit placed, tandem placements kept, ALT clip kept on the ALT, byte-identical across -t, --fast and --dedup-reads)"
