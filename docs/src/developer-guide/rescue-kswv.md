# The 11-op rescue cell

The rescue kernels in `src/kswv.cpp` run many Smith-Waterman jobs side by side, one per SIMD lane, 8- or 16-bit. Every SIMD body can build its cell two ways, and `BWA3_RESCUE_FSCAN` picks one: the original cell, or the cell rebuilt around the score before the in-row gap (the "FScan" cell, default on). The NEON 8-bit kernel is where it was first written and where the numbers in the name come from; the same cell is in the NEON 16-bit body (`kswv_neon_16`) and the AVX2 and AVX-512BW bodies (`kswv256_u8`, `kswv256_16`, `kswv512_u8`, `kswv512_16`).

## The cell

The original cell takes H as the maximum of the diagonal score M, the in-row gap E and the vertical gap F, then opens both next gaps from H. The fused cell first forms G = max(M, F), the score before the in-row gap, takes H = max(G, E), and opens both next gaps from one saturating subtract t = G - (o + e), extending them from E and F as before.

- Opening from G instead of H drops two transitions: an in-row gap run followed directly by an in-row gap open, which is never better than extending the run, and an in-row gap run followed directly by a vertical one, which always has an equal-scoring twin with the two runs swapped. So every H the kernel compares, and with it every score and position it reports, is unchanged. The argument is written out at `KSWV_NEON_U8_CELL_PAIR_FS` (16-bit: `KSWV_NEON_16_CELL_FS`). Gates: `kswv u8 rescue: BWA3_RESCUE_FSCAN off == on in every u8 body, and FSCAN matches scalar`, `kswv u16 rescue: BWA3_RESCUE_FSCAN off == on in every u16 body, and FSCAN matches scalar`, `kswv --meth freed cell: BWA3_RESCUE_FSCAN off == on, and FSCAN matches scalar`.
- One subtract serves both gaps only when the insertion and deletion open-plus-extend sums are equal; the opens themselves may differ. The sum must fit the lane (255 for 8-bit, `INT16_MAX` for 16-bit), and the argument needs every open and extend to be non-negative. `fscan_scoring_ok` checks all of it for every dispatcher, and any other scoring takes the original body. bwa's defaults and presets all qualify. Gate: `kswv rescue: the FSCAN gate keeps asymmetric gaps exact`.
- The cell needs no reference-padding mask: a pad row's diagonal term cannot exceed the score it comes from, as the note above `kswv_neon_u8` shows, so only the query half of the boundary mask remains, from the column where the lanes' query ends begin (`compute_jsplit`).

In the NEON 8-bit kernel the cell is 11 vector operations instead of 13, and the next row's cell no longer waits on this row's in-row gap chain, which is what lets the two-row sweep overlap them. The other bodies save one max or subtract per cell and the pad mask.

The banded rescue DP uses the same fused cell ([banded rescue DP](rescue-banding.md#the-kernels)); its band argument has to cover the dropped transition separately, and does.

## The other toggles

Three older switches each choose between two forms of a kernel step, with a unit test per switch against the scalar reference. All three are NEON switches; `BWA3_RESCUE_USQADD` also reaches the AVX2 8-bit FScan body:

- `BWA3_RESCUE_USQADD`: in the 8-bit kernel, one saturating add of a signed score delta per cell instead of the biased add and de-biasing subtract. NEON has an unsigned-plus-signed saturating add (USQADD). x86 has none, so the AVX2 FScan body instead keeps H, E and F as H - 128 in signed bytes: that map is an order-preserving bijection of 0..255 onto -128..127 that takes unsigned saturation at 0 to signed saturation at -128, so every step of the cell has an exact signed twin and the score becomes one signed saturating add. The row maximum goes back to the unsigned domain once per row. It needs every gap constant to fit a positive signed byte: the open-plus-extend sum at most 127, which bounds each extend too, since the FScan gate already requires non-negative costs. Above that the biased FScan body runs. Gate: `kswv u8 rescue: BWA3_RESCUE_USQADD off == on, and biased body matches scalar`.
- `BWA3_RESCUE_ROWPAIR`: sweep two target rows per pass over the query, sharing the column loads and the vertical carry. Gates: `kswv u8 rescue: BWA3_RESCUE_ROWPAIR off == on, and one-row matches scalar`, `kswv u16 rescue: ROWPAIR/LAZYQE configurations agree, and one-row matches scalar`.
- `BWA3_RESCUE_LAZYQE`: in the two-row sweep, find each row's query end from the stored row after the sweep instead of tracking it inline; both forms take the first column holding the row maximum, so the result is identical. Gate: `kswv u8 rescue: BWA3_RESCUE_LAZYQE off == on in the two-row sweep, and inline matches scalar`. The AVX2 16-bit body always finds its query end this way (no switch): it checkpoints the running row maximum once per block of columns and scans only the blocks that can hold it. Gates: `kswv::getScores16 matches scalar ksw_align2 on 10k random + curated edge pairs`, `kswv::getScores16 matches scalar on mixed-query-length batches`, `kswv u16 rescue: BWA3_RESCUE_FSCAN off == on in every u16 body, and FSCAN matches scalar`.

These toggles are read on every kernel call rather than once, so a unit test can flip them in-process; the cost is one lookup per lane group of a batch, never per cell.

## Coverage by tier

The unit tests run on the NEON and AVX2 CI rows. The AVX-512BW bodies have no unit-test row, since a CI runner may lack the instructions; they are checked end to end by `rescue_prune_identity.sh` under `BWAMEM3_FORCE_TIER`, on the canonical row when its runner has AVX-512BW, and reported skipped when it does not. Gates: `Rescue identity under each forced x86 tier (generated PE fixture)`, `rescue_prune_identity.sh`.
