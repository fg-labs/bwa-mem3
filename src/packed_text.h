#ifndef BWA_PACKED_TEXT_H
#define BWA_PACKED_TEXT_H

#include <cstddef>
#include <cstdint>
#include <string>

// Read-only mmap wrapper over a `.pac`-format packed text.
// Text layout matches bwa-mem2's _set_pac/_get_pac macros (see bntseq.cpp):
// base at position l occupies bits ((~l & 3) << 1)+1 .. ((~l & 3) << 1)
// of byte pac[l >> 2]. MSB-first within a byte.
class PackedText {
public:
    // Opens `pac_path` read-only, mmaps its contents. `n_bases` is the
    // logical text length in bases. The file must contain at least
    // ceil(n_bases / 4) bytes.
    PackedText(const std::string& pac_path, int64_t n_bases);
    ~PackedText();

    PackedText(const PackedText&) = delete;
    PackedText& operator=(const PackedText&) = delete;

    int64_t       length() const { return n_; }
    const uint8_t* data() const  { return data_; }

    // Return the 2-bit base at position `pos`. Out-of-range positions return 0.
    // Inlined so the hot preludes (emit_0123, compute_counts) vectorise.
    uint8_t get_base(int64_t pos) const {
        if (pos < 0 || pos >= n_) return 0;
        int shift = (int)((~pos & 3) << 1);
        return (uint8_t)((data_[pos >> 2] >> shift) & 3);
    }

    // Return up to `k` bases (k <= 32) starting at `pos`, packed MSB-first;
    // bases past end-of-text are zero-padded (sentinel-small).
    uint64_t get_kmer(int64_t pos, int k) const;

private:
    int     fd_ = -1;
    size_t  map_bytes_ = 0;
    const uint8_t* data_ = nullptr;
    int64_t n_ = 0;
};

#endif
