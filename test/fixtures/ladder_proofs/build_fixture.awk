# Deterministically construct the reference and the single-end reads for
# test/ladder_rungs_test.sh and test/ungapped_hit_parity_test.sh. All bases
# are sliced from committed phix.fa, and every substitution is a fixed
# rotation (A>C>G>T>A), so the fixture is byte-identical across awk
# implementations -- no PRNG.
#
# Run as:
#   mawk -v MODE=ref               -f build_fixture.awk test/fixtures/phix.fa > ref.fa
#   mawk -v MODE=read -v CASE=<c>  -f build_fixture.awk test/fixtures/phix.fa > read.fq
#
# The reference is phix[0, 2650) (0-based, half-open) under one name. Every
# read below is anchored at ref[1000, ...).
#
# CASE=ladder -- the retry-ladder rung count (upstream MAX_BAND_TRY 2).
#   Reference layout used by the read:
#     S   = phix[1000,1210)  210 bp, the read's only exact seed
#     Q2a = phix[1210,1350)  140 bp
#     Q2b = phix[1350,1480)  130 bp
#     Q3  = phix[1480,1650)  170 bp
#   Read (900 bp): S | I1 | mut(Q2a) | I2 | mut(Q2b) | I3 | mut(Q3), where
#     I1 = phix[3000,3080) 80 bp, I2 = phix[3100,3170) 70 bp,
#     I3 = phix[3200,3300) 100 bp are insertions relative to the reference
#     (their bases lie outside the reference), and mut() substitutes every
#     18th base so the segment holds no exact 19-mer and seeds nothing.
#   Right extension from S at default scoring: each Q segment nets +19 over
#   its insertion (140-7*5-86, 130-7*5-76, 170-9*5-106) at diagonal offsets
#   80, 150 and 250. The ladder at w=100 finds Q2a (max_off 80 >= 75, score
#   changed), at 2w=200 finds Q2b (max_off 150 >= 150, score changed).
#   Upstream stops there: 209M80I141M70I125M275S, AS 248 (bwa 0.7.19). A
#   third rung at 400 reaches Q3 (offset 250): 209M80I141M70I130M100I170M,
#   AS 267.
#
# CASE=adaptive -- the --adaptive-band narrowing ladder's rung count
#   (ADAPTIVE_BAND_TRY 4: band_start 20 << i on the 16-bit tier).
#   Read (554 bp): phix[1000,1100) seed | I 16 | mut(phix[1100,1220)) | I 16 |
#     mut(phix[1220,1340)) | I 32 | mut(phix[1340,1490)), insertions from
#     phix[3400,...). The right extension's gains sit at diagonal offsets 16, 32
#     and 64, each just past 3/4 of the rung (20, 40, 80) that first sees it, so
#     a 2-rung narrowing ladder stops at 40 and a 4-rung one reaches the third.
#
# CASE=ladder_rc -- the same read reverse-complemented. It maps to the reverse
#   strand, so the seed sits at the 3' end of the query and the staircase is
#   extended to the LEFT of it: the left-extension ladder loops.
#
# CASE=L0tie, O20tie, d3zdrop -- the ungapped HIT path vs the extension
#   kernel. Each read is phix[1000,1150) (150 bp) with substitutions in its
#   first bases; the exact seed is the rest of the read and the LEFT extension
#   walks the reversed prefix away from it.
#   L0tie  : read[4]         -> reversed prefix [X][m m m m]. Under -L 0 the
#            diagonal returns to h0 without exceeding it; the kernels keep the
#            earlier record (qle 0), so bwa clips: 5S145M.
#   O20tie : read[0,1,6]     -> reversed prefix [X][m m m m][X][X] under -O 20
#            (x_threshold 4): tie at step 5, then -8; branch A at the default
#            -L 5, kernels' qle 0: 7S143M.
#   d3zdrop: read[10]        -> reversed prefix [X][m x10] under -d 3: the
#            kernels z-drop on the mismatch (drop 4 > 3): 11S139M, AS 139.

/^>/ { next }
     { phix = phix $0 }

function comp(b) {
    if (b == "A") return "T"
    if (b == "C") return "G"
    if (b == "G") return "C"
    if (b == "T") return "A"
    return "N"
}

# Reverse complement of s.
function revcomp(s,    out, i) {
    out = ""
    for (i = length(s); i >= 1; i--) out = out comp(substr(s, i, 1))
    return out
}

function rot(b) {
    if (b == "A") return "C"
    if (b == "C") return "G"
    if (b == "G") return "T"
    return "A"
}

# Substitute base at 1-based position p of s.
function subst(s, p) {
    return substr(s, 1, p - 1) rot(substr(s, p, 1)) substr(s, p + 1)
}

# Substitute every `every`-th base of s (positions every, 2*every, ...).
function mut(s, every,    out, i) {
    out = s
    for (i = every; i <= length(s); i += every) out = subst(out, i)
    return out
}

function emit_fa(name, s,    i) {
    printf ">%s\n", name
    for (i = 1; i <= length(s); i += 70) printf "%s\n", substr(s, i, 70)
}

function emit_fq(name, s,    q, i) {
    q = ""
    for (i = 1; i <= length(s); i++) q = q "I"
    printf "@%s\n%s\n+\n%s\n", name, s, q
}

END {
    if (length(phix) < 3300) {
        print "ERROR: phix body too short (" length(phix) " < 3300 bp)" > "/dev/stderr"
        exit 1
    }
    if (MODE == "ref") {
        emit_fa("ladder_proofs_ref", substr(phix, 1, 2650))
        exit 0
    }
    if (MODE != "read") {
        print "ERROR: MODE must be ref or read" > "/dev/stderr"
        exit 1
    }
    if (CASE == "ladder" || CASE == "ladder_rc") {
        s   = substr(phix, 1001, 210)
        q2a = substr(phix, 1211, 140)
        q2b = substr(phix, 1351, 130)
        q3  = substr(phix, 1481, 170)
        i1  = substr(phix, 3001, 80)
        i2  = substr(phix, 3101, 70)
        i3  = substr(phix, 3201, 100)
        read = s i1 mut(q2a, 18) i2 mut(q2b, 18) i3 mut(q3, 18)
        if (length(read) != 900) {
            print "ERROR: ladder read is " length(read) " bp, expected 900" > "/dev/stderr"
            exit 1
        }
        if (CASE == "ladder_rc") read = revcomp(read)
        emit_fq(CASE, read)
        exit 0
    }
    if (CASE == "adaptive") {
        read = substr(phix, 1001, 100) substr(phix, 3401, 16) mut(substr(phix, 1101, 120), 18) \
               substr(phix, 3421, 16) mut(substr(phix, 1221, 120), 18) \
               substr(phix, 3441, 32) mut(substr(phix, 1341, 150), 18)
        if (length(read) != 554) {
            print "ERROR: adaptive read is " length(read) " bp, expected 554" > "/dev/stderr"
            exit 1
        }
        emit_fq("adaptive", read)
        exit 0
    }
    read = substr(phix, 1001, 150)
    if (CASE == "L0tie") {
        read = subst(read, 5)
    } else if (CASE == "O20tie") {
        read = subst(read, 1); read = subst(read, 2); read = subst(read, 7)
    } else if (CASE == "d3zdrop") {
        read = subst(read, 11)
    } else {
        print "ERROR: unknown CASE " CASE > "/dev/stderr"
        exit 1
    }
    emit_fq(CASE, read)
}
