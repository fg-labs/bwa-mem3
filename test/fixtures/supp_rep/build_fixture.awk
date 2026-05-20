# Deterministically construct an engineered reference + SE reads that
# trigger --supp-rep-hard-cap. All bases are sliced from committed
# phix.fa so the generated fixture is byte-identical across awk
# implementations (mawk version drift, gawk, BSD awk) — no PRNG.
#
# Run as:
#   mawk -v MODE=ref   -f build_fixture.awk test/fixtures/phix.fa > ref.fa
#   mawk -v MODE=reads -f build_fixture.awk test/fixtures/phix.fa > reads.fq
#
# Phix slice allocations (positions 0-based, half-open; phix.fa is 5386 bp):
#   [   0,  400)  UNIQUE_A   (400 bp, 1 copy in ref)
#   [ 400,  600)  UNIQUE_C   (200 bp, 1 copy in ref; disambiguates supp chain)
#   [ 600,  660)  MOTIF      ( 60 bp, 8 copies in ref; high SA-count seed)
#   [ 700, 4100)  filler bank — see allocations below
#
# Filler slot allocations (all from phix, each used exactly once):
#   F_pre        (100 bp, phix[ 700, 800))   — between ref start and UNIQUE_A
#   F_post       (200 bp, phix[ 800,1000))   — between UNIQUE_A and MOTIF copy 0
#   F_inter[k]   (500 bp, phix[1000+500*k, 1500+500*k)) — between adjacent
#                                                          MOTIF copies (6 slots)
#   F_tail       (100 bp, phix[4000,4100))   — after the last MOTIF copy
#
# Ref layout (each | separates one concatenated piece):
#   F_pre | UA | F_post | M | F_i0 | M | F_i1 | M | F_i2 | M | F_i3 | M | UC
#                                                            \-- DISAMBIG (5th M)
#         | M | F_i4 | M | F_i5 | M | F_tail
#
# Eight MOTIF copies total; the fifth (0-indexed DISAMBIG=4) is followed
# immediately by UNIQUE_C, so the 260 bp string MOTIF+UNIQUE_C appears
# uniquely once in the ref while MOTIF alone appears eight times. The
# 500 bp inter-repeat gap is large enough that bwa-mem2 emits each MOTIF
# copy as a separate chain seed (SA-count 8 → n_hits=8); too small a gap
# can let chain merging collapse the eight ref hits into a single chain
# without propagating the per-seed n_hits=8 to chain_n_hits.
#
# Read structure: UNIQUE_A + MOTIF + UNIQUE_C = 660 bp.
#   primary chain → UNIQUE_A   (chain_n_hits = 1, MAPQ = 60)
#   supp chain    → MOTIF + UNIQUE_C at DISAMBIG  (chain_n_hits = 8 from the
#                                                  MOTIF seed; natural MAPQ
#                                                  high because UNIQUE_C
#                                                  disambiguates the chain)
# With cap = 0 the supp keeps its natural MAPQ. With cap >= 2 (and <= 8) the
# supp goes to MAPQ = 0.

# Accumulate phix bases (concatenate every non-header line).
/^>/ { next }
     { phix = phix $0 }

END {
    if (length(phix) < 4100) {
        print "ERROR: phix body too short (" length(phix) " < 4100 bp)" > "/dev/stderr"
        exit 1
    }

    # 1-based substr offsets. Phix position k (0-based) is substr(phix, k+1).
    unique_a = substr(phix, 1,   400)        # phix[0,   400)
    unique_c = substr(phix, 401, 200)        # phix[400, 600)
    motif    = substr(phix, 601,  60)        # phix[600, 660)

    f_pre    = substr(phix, 701, 100)        # phix[ 700,  800)
    f_post   = substr(phix, 801, 200)        # phix[ 800, 1000)
    f_tail   = substr(phix, 4001, 100)       # phix[4000, 4100)
    for (k = 0; k < 6; k++) {
        f_inter[k] = substr(phix, 1001 + k * 500, 500)  # phix[1000 + 500k, 1500 + 500k)
    }

    N_COPIES      = 8
    DISAMBIG_COPY = 4   # 0-indexed copy after which UNIQUE_C is inserted

    ref = f_pre unique_a f_post
    inter_idx = 0
    for (i = 0; i < N_COPIES; i++) {
        ref = ref motif
        if (i == DISAMBIG_COPY)    ref = ref unique_c
        else if (i < N_COPIES - 1) { ref = ref f_inter[inter_idx]; inter_idx++ }
    }
    ref = ref f_tail
    L = length(ref)

    if (MODE == "ref") {
        print ">supp_rep_engineered"
        for (i = 1; i <= L; i += 70) print substr(ref, i, 70)
        exit
    }
    if (MODE == "reads") {
        # 30 SE reads, identical body so the cap effect is uniform; distinct
        # read names so SAM record counts are stable.
        read = unique_a motif unique_c
        qual = ""
        for (q = 0; q < length(read); q++) qual = qual "I"
        for (n = 0; n < 30; n++) printf "@r%02d\n%s\n+\n%s\n", n, read, qual
        exit
    }
    print "ERROR: MODE must be ref or reads" > "/dev/stderr"
    exit 1
}
