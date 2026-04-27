/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************************************/

#include "fm_index_writer.h"
#include "FMI_search.h"   // CP_OCC, CP_BLOCK_SIZE, CP_SHIFT, CP_MASK, DUMMY_CHAR
#include "io_utils.h"     // pwrite_all
#include "macro.h"        // SA_COMPRESSION, SA_COMPX, SA_COMPX_MASK
#include "utils.h"        // err_fatal

#include <algorithm>
#include <array>
#include <cerrno>
#include <string>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// A small ring buffer that accumulates writes to one output region. When
// `fill` reaches `cap`, the whole ring is pwritten at `cursor` and `cursor`
// advances. This amortises the syscall cost across many small appends.
struct Ring {
    int      fd;
    uint8_t* buf;
    size_t   cap;      // bytes
    size_t   fill;     // bytes used
    off_t    cursor;   // next file offset to write
    const char* what;  // for error messages
};

static inline void ring_flush(Ring& r) {
    if (r.fill == 0) return;
    pwrite_all(r.fd, r.buf, r.fill, r.cursor, r.what);
    r.cursor += (off_t)r.fill;
    r.fill    = 0;
}

static inline void ring_append(Ring& r, const void* src, size_t n) {
    if (r.fill + n <= r.cap) {
        memcpy(r.buf + r.fill, src, n);
        r.fill += n;
        return;
    }
    ring_flush(r);
    if (n >= r.cap) {
        pwrite_all(r.fd, src, n, r.cursor, r.what);
        r.cursor += (off_t)n;
        return;
    }
    memcpy(r.buf + r.fill, src, n);
    r.fill += n;
}

// Build a CP_OCC entry from a 64-byte BWT window (may be padded with
// DUMMY_CHAR for a partial final block) and a snapshot of cumulative counts
// taken at the start of the block. Bit 63-j of one_hot_bwt_str[c] is 1 iff
// bwt_block[j] == c. Equivalent to the older "<<= 1; += 1" formulation but
// without the 256-cycle serial shift chain.
static inline CP_OCC make_cp_occ(const uint8_t bwt_block[CP_BLOCK_SIZE],
                                 const int64_t cp_count_snap[4]) {
    CP_OCC cpo;
    cpo.cp_count[0]        = cp_count_snap[0];
    cpo.cp_count[1]        = cp_count_snap[1];
    cpo.cp_count[2]        = cp_count_snap[2];
    cpo.cp_count[3]        = cp_count_snap[3];
    cpo.one_hot_bwt_str[0] = 0;
    cpo.one_hot_bwt_str[1] = 0;
    cpo.one_hot_bwt_str[2] = 0;
    cpo.one_hot_bwt_str[3] = 0;
    for (int j = 0; j < CP_BLOCK_SIZE; ++j) {
        uint8_t c = bwt_block[j];
        if (c < 4) cpo.one_hot_bwt_str[c] |= (uint64_t)1 << (63 - j);
    }
    return cpo;
}

// Fetch BWT byte for SA row i in bwa's alphabet (A=0..T=3, sentinel=4),
// given SA[i] = pos and the libsais-alphabet input buffer (0=$, 1..4=ACGT).
// pos == 0 means the full-text suffix; its predecessor is the GSA
// terminator at buf[N] (value 0, which maps to bwa's sentinel).
static inline uint8_t bwt_byte(int64_t pos, int64_t N, const uint8_t* buf) {
    int64_t src = (pos == 0) ? N : pos - 1;
    uint8_t c   = buf[src];
    return (c == 0) ? 4 : (uint8_t)(c - 1);
}

// Pass 1: count BWT bytes in stripe [s, e) and locate the sentinel row.
template <typename SaIdx>
static void count_stripe(int64_t s, int64_t e,
                         const uint8_t* buf, const SaIdx* sa, int64_t N,
                         int64_t out_count[5],
                         int64_t& out_local_sentinel)
{
    int64_t cnt[5] = {0, 0, 0, 0, 0};
    int64_t local_sentinel = -1;
    for (int64_t i = s; i < e; ++i) {
        uint8_t b = bwt_byte((int64_t)sa[i], N, buf);
        ++cnt[b];
        if (b == 4) local_sentinel = i;
    }
    for (int k = 0; k < 5; ++k) out_count[k] = cnt[k];
    out_local_sentinel = local_sentinel;
}

// Pass 2: emit cp_occ blocks and SA samples for stripe [s, e). `s` is
// always CP_BLOCK_SIZE-aligned (by stripe_start construction); `e` equals
// ref_seq_len only in the last stripe and may leave a partial block.
// `running_init` is the cumulative base count at row s from pass 1's
// prefix sum. Ring buffers pwrite into disjoint file ranges.
//
// Structure: outer loop over full 64-row blocks, inner branchless 64-row
// inner loop + a single cp_occ emit per block. The per-block snap of
// running counts + per-8-row sample emit are hoisted out of the inner
// loop so it has no per-row branches. Partial final block handled after.
template <typename SaIdx>
static void emit_stripe(int64_t s, int64_t e,
                        const uint8_t* buf, const SaIdx* sa, int64_t N,
                        int fd,
                        off_t off_cp, off_t off_msb, off_t off_lsw,
                        const int64_t running_init[4])
{
    constexpr size_t kSampleRingSamples = 4096;
    constexpr size_t kCpOccRingEntries  = 64;
    std::vector<uint8_t> cp_occ_buf(kCpOccRingEntries * sizeof(CP_OCC));
    std::vector<uint8_t> sa_ms_buf(kSampleRingSamples * sizeof(int8_t));
    std::vector<uint8_t> sa_ls_buf(kSampleRingSamples * sizeof(uint32_t));
    Ring cp_occ_ring = { fd, cp_occ_buf.data(), cp_occ_buf.size(), 0, off_cp,  "cp_occ" };
    Ring sa_ms_ring  = { fd, sa_ms_buf.data(),  sa_ms_buf.size(),  0, off_msb, "sa_ms_byte" };
    Ring sa_ls_ring  = { fd, sa_ls_buf.data(),  sa_ls_buf.size(),  0, off_lsw, "sa_ls_word" };

    // running[5]: bases A/C/G/T occurrence-count accumulators ([0..3]) plus
    // a trailing sentinel-bucket ([4]) that absorbs the b==4 increment from
    // the inner loop's `running[b]++`. The sentinel slot is never read --
    // cp_count_snap only copies [0..3] -- but keeping it lets the hot inner
    // loop stay branchless (no `if (b < 4)` guard before the increment).
    int64_t running[5] = { running_init[0], running_init[1],
                           running_init[2], running_init[3], 0 };

    // Full blocks: s is CP_BLOCK_SIZE-aligned, so blk_start is too, and
    // since CP_BLOCK_SIZE (64) is a multiple of SA_COMPX (8) the global
    // (i & SA_COMPX_MASK) check reduces to a block-local (j & ...) check.
    int64_t full_end = s + ((e - s) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    uint8_t bwt_block[CP_BLOCK_SIZE];

    for (int64_t blk_start = s; blk_start < full_end; blk_start += CP_BLOCK_SIZE) {
        int64_t cp_count_snap[4] = { running[0], running[1], running[2], running[3] };
        for (int j = 0; j < CP_BLOCK_SIZE; ++j) {
            int64_t pos = (int64_t)sa[blk_start + j];
            uint8_t b   = bwt_byte(pos, N, buf);
            bwt_block[j] = b;
            running[b]++;
            if ((j & SA_COMPX_MASK) == 0) {
                int8_t   ms = (int8_t)(((uint64_t)pos >> 32) & 0xFF);
                uint32_t ls = (uint32_t)(pos & 0xFFFFFFFFULL);
                ring_append(sa_ms_ring, &ms, sizeof(ms));
                ring_append(sa_ls_ring, &ls, sizeof(ls));
            }
        }
        CP_OCC cpo = make_cp_occ(bwt_block, cp_count_snap);
        ring_append(cp_occ_ring, &cpo, sizeof(CP_OCC));
    }

    // Partial final block (last stripe only, when ref_seq_len isn't a
    // multiple of CP_BLOCK_SIZE).
    if (full_end < e) {
        int64_t cp_count_snap[4] = { running[0], running[1], running[2], running[3] };
        int remaining = (int)(e - full_end);
        for (int j = 0; j < remaining; ++j) {
            int64_t pos = (int64_t)sa[full_end + j];
            uint8_t b   = bwt_byte(pos, N, buf);
            bwt_block[j] = b;
            running[b]++;
            if ((j & SA_COMPX_MASK) == 0) {
                int8_t   ms = (int8_t)(((uint64_t)pos >> 32) & 0xFF);
                uint32_t ls = (uint32_t)(pos & 0xFFFFFFFFULL);
                ring_append(sa_ms_ring, &ms, sizeof(ms));
                ring_append(sa_ls_ring, &ls, sizeof(ls));
            }
        }
        for (int k = remaining; k < CP_BLOCK_SIZE; ++k)
            bwt_block[k] = DUMMY_CHAR;
        CP_OCC cpo = make_cp_occ(bwt_block, cp_count_snap);
        ring_append(cp_occ_ring, &cpo, sizeof(CP_OCC));
    }

    ring_flush(cp_occ_ring);
    ring_flush(sa_ms_ring);
    ring_flush(sa_ls_ring);
}

// Type-specialised driver: pass 1 counting + pass 2 emit over the supplied
// SA element type. The public entry point erases the type via sa_is_64bit.
// Returns nothing: every failure path inside calls err_fatal and terminates
// the process, so a recoverable-error int rc would be a misleading contract.
template <typename SaIdx>
static void write_fm_index_streaming_typed(int fd,
                                           const uint8_t* buf,
                                           const SaIdx* sa,
                                           int64_t pac_len,
                                           int64_t ref_seq_len,
                                           off_t off_cp_occ,
                                           off_t off_ms_byte,
                                           off_t off_ls_word,
                                           off_t off_sent,
                                           const int64_t count[5],
                                           int64_t* out_sentinel_index,
                                           int T)
{
    std::vector<int64_t> stripe_start(T + 1, 0);
    for (int t = 0; t <= T; ++t) {
        int64_t pos = (ref_seq_len * (int64_t)t) / T;
        pos -= pos % CP_BLOCK_SIZE;
        stripe_start[t] = pos;
    }
    stripe_start[T] = ref_seq_len;

    // count_stripe / emit_stripe rely on stripe_start being monotonically
    // non-decreasing so adjacent stripes index disjoint half-open ranges.
    // The integer-arithmetic + CP_BLOCK_SIZE-rounding step above preserves
    // monotonicity, but stripe_start[T] is overwritten with the unrounded
    // ref_seq_len; assert the invariant rather than letting a future
    // regression silently produce a corrupt FM index.
    for (int t = 0; t < T; ++t) {
        if (stripe_start[t] > stripe_start[t + 1])
            err_fatal("write_fm_index_streaming",
                      "stripe_start non-monotonic: stripe_start[%d]=%lld > stripe_start[%d]=%lld "
                      "(ref_seq_len=%lld, CP_BLOCK_SIZE=%lld, T=%d)",
                      t, (long long)stripe_start[t],
                      t + 1, (long long)stripe_start[t + 1],
                      (long long)ref_seq_len, (long long)CP_BLOCK_SIZE, T);
    }

    // Pass 1: per-stripe BWT-base histograms + local sentinel.
    std::vector<std::array<int64_t, 5>> stripe_counts(T);
    std::vector<int64_t> stripe_sentinel(T, -1);
#ifdef _OPENMP
    #pragma omp parallel for num_threads(T) schedule(static)
#endif
    for (int t = 0; t < T; ++t) {
        int64_t local_sentinel = -1;
        count_stripe<SaIdx>(stripe_start[t], stripe_start[t + 1],
                            buf, sa, pac_len,
                            stripe_counts[t].data(), local_sentinel);
        stripe_sentinel[t] = local_sentinel;
    }

    // Reduce sentinel + prefix-sum running counts at stripe starts.
    int64_t sentinel_index = -1;
    for (int t = 0; t < T; ++t) {
        if (stripe_sentinel[t] < 0) continue;
        if (sentinel_index >= 0)
            err_fatal("write_fm_index_streaming",
                      "multiple sentinels in BWT (at %lld and %lld)",
                      (long long)sentinel_index, (long long)stripe_sentinel[t]);
        sentinel_index = stripe_sentinel[t];
    }
    if (sentinel_index < 0)
        err_fatal("write_fm_index_streaming", "BWT has no sentinel");

    // Belt-and-braces: count_stripe records only the last sentinel per
    // stripe, so two sentinels within the same stripe would be silently
    // collapsed. The cross-stripe reduce above only catches duplicates
    // that fall in different stripes. Sum the per-stripe sentinel
    // histograms to enforce the global "exactly one sentinel" invariant.
    int64_t total_sentinels = 0;
    for (int t = 0; t < T; ++t) total_sentinels += stripe_counts[t][4];
    if (total_sentinels != 1)
        err_fatal("write_fm_index_streaming",
                  "BWT has %lld sentinels, expected 1",
                  (long long)total_sentinels);

    std::vector<std::array<int64_t, 4>> stripe_running_start(T);
    {
        int64_t acc[4] = {0, 0, 0, 0};
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < 4; ++k) stripe_running_start[t][k] = acc[k];
            for (int k = 0; k < 4; ++k) acc[k] += stripe_counts[t][k];
        }
    }

    // Pass 2: per-stripe emit to disjoint file ranges.
#ifdef _OPENMP
    #pragma omp parallel for num_threads(T) schedule(static)
#endif
    for (int t = 0; t < T; ++t) {
        int64_t s = stripe_start[t];
        int64_t e = stripe_start[t + 1];
        if (s >= e) continue;

        int64_t first_block_idx  = s >> CP_SHIFT;
        int64_t first_sample_idx = s >> SA_COMPX;
        off_t off_cp  = off_cp_occ  + (off_t)first_block_idx  * (off_t)sizeof(CP_OCC);
        off_t off_msb = off_ms_byte + (off_t)first_sample_idx * (off_t)sizeof(int8_t);
        off_t off_lsw = off_ls_word + (off_t)first_sample_idx * (off_t)sizeof(uint32_t);

        emit_stripe<SaIdx>(s, e, buf, sa, pac_len,
                           fd, off_cp, off_msb, off_lsw,
                           stripe_running_start[t].data());
    }

    pwrite_all(fd, &ref_seq_len,     sizeof(int64_t),     0,               "ref_seq_len");
    pwrite_all(fd, count,            5 * sizeof(int64_t), sizeof(int64_t), "count[5]");
    pwrite_all(fd, &sentinel_index,  sizeof(int64_t),     off_sent,        "sentinel_index");

    *out_sentinel_index = sentinel_index;
}

void write_fm_index_streaming(const char* out_path,
                              const uint8_t* buf,
                              const void* sa,
                              bool sa_is_64bit,
                              int64_t pac_len,
                              const int64_t count[5],
                              int64_t* out_sentinel_index,
                              int num_threads)
{
    const int64_t ref_seq_len = pac_len + 1;
    // sa_ms_byte is stored as int8_t (top byte of a 40-bit SA value), so
    // SA values >= 2^39 would silently sign-extend on load and corrupt the
    // index. Reject up front rather than producing a broken .bwt.2bit.64.
    if (ref_seq_len > ((int64_t)1 << 39))
        err_fatal(__func__,
                  "ref_seq_len=%lld exceeds on-disk SA-sample range (~5.5e11)",
                  (long long)ref_seq_len);

    const off_t HDR_BYTES       = (off_t)(sizeof(int64_t) + 5 * sizeof(int64_t));
    const int64_t cp_occ_size   = (ref_seq_len >> CP_SHIFT) + 1;
    const int64_t sa_sample_cnt = (ref_seq_len >> SA_COMPX) + 1;

    const off_t off_cp_occ  = HDR_BYTES;
    const off_t off_ms_byte = off_cp_occ  + (off_t)cp_occ_size   * (off_t)sizeof(CP_OCC);
    const off_t off_ls_word = off_ms_byte + (off_t)sa_sample_cnt * (off_t)sizeof(int8_t);
    const off_t off_sent    = off_ls_word + (off_t)sa_sample_cnt * (off_t)sizeof(uint32_t);
    const off_t total_bytes = off_sent    + (off_t)sizeof(int64_t);

    // Write to a sibling .tmp file and rename on success so a failed
    // build (err_fatal anywhere downstream) never leaves a partial /
    // zero-filled index sitting at the canonical path. A leftover .tmp
    // from a previous failed run is overwritten by O_TRUNC on the next
    // attempt, so it self-heals without explicit cleanup.
    const std::string tmp_path = std::string(out_path) + ".tmp";
    int fd = open(tmp_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0)
        err_fatal(__func__, "open('%s'): %s", tmp_path.c_str(), strerror(errno));
    if (ftruncate(fd, total_bytes) != 0)
        err_fatal(__func__, "ftruncate('%s', %lld): %s",
                  tmp_path.c_str(), (long long)total_bytes, strerror(errno));

    // For tiny references (full_blocks < 2 * T), per-stripe ring buffers
    // and OMP pool spin-up don't amortise, so we silently scale T down,
    // collapsing to 1 when full_blocks < 2. Callers asking for many
    // threads on a sub-CP_BLOCK_SIZE input may therefore see only one
    // stripe in the streaming writer's logs — that is intentional.
    int T = num_threads > 0 ? num_threads : 1;
    if (T > 1) {
        int64_t full_blocks = ref_seq_len / CP_BLOCK_SIZE;
        if (full_blocks < T * 2) T = (int)std::max<int64_t>(1, full_blocks / 2);
    }

    if (sa_is_64bit) {
        write_fm_index_streaming_typed<int64_t>(
            fd, buf, (const int64_t*)sa, pac_len, ref_seq_len,
            off_cp_occ, off_ms_byte, off_ls_word, off_sent,
            count, out_sentinel_index, T);
    } else {
        write_fm_index_streaming_typed<int32_t>(
            fd, buf, (const int32_t*)sa, pac_len, ref_seq_len,
            off_cp_occ, off_ms_byte, off_ls_word, off_sent,
            count, out_sentinel_index, T);
    }

    if (close(fd) != 0)
        err_fatal(__func__, "close('%s'): %s", tmp_path.c_str(), strerror(errno));
    // Atomic publish: rename(2) is atomic on the same filesystem, so a
    // concurrent reader sees either the previous good index or the new
    // one, never a partial.
    if (rename(tmp_path.c_str(), out_path) != 0)
        err_fatal(__func__, "rename('%s' -> '%s'): %s",
                  tmp_path.c_str(), out_path, strerror(errno));
}
