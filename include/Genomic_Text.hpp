
#ifndef CAPS_SA_GENOMIC_TEXT_HPP
#define CAPS_SA_GENOMIC_TEXT_HPP



#include <cstdint>
#include <cstddef>
#include <vector>
#include <cassert>


namespace CaPS_SA
{

// Class to to represent genomic texts in a 2-bit-packed manner. The bit-packed
// representation is considered to run right-to-left.
class Genomic_Text
{
    const std::size_t n_;   // Length of the input text.

    const std::size_t pack_sz;  // Size of the bit-packed text.
    std::vector<uint8_t> B; // The bit-packed text.

    // Returns the 124-nucleobase block (31 bytes) from onward the `i`'th
    // nucleobase, in 256-bits little-endian. No guarantees are provided for
    // the highest byte.
    __m256i load(std::size_t i) const;

public:

    // Constructs a 2-bit-packed representation of the genomic text `T` of
    // length `n`.
    Genomic_Text(const char* T, std::size_t n);

    Genomic_Text(const Genomic_Text&)  = delete;
    Genomic_Text(Genomic_Text&&) = delete;
    Genomic_Text& operator=(const Genomic_Text&) = delete;
    Genomic_Text& operator=(Genomic_Text&&) = delete;

    // Returns the code of the nucleobase at index `idx` of the original text.
    uint8_t base_at(std::size_t idx) const;
};


inline uint8_t Genomic_Text::base_at(const std::size_t idx) const
{
    assert(idx < n_);

    const auto byte_idx = idx >> 2;
    const auto bit_idx = (idx & 0b11) * 2;

    return  (B[byte_idx] & (0b11 << bit_idx)) >> bit_idx;
}


inline __m256i Genomic_Text::load(const std::size_t i) const
{
    assert(i + 124 <= n_);

    const auto base = i / 4;    // Base word's index.
    const auto blk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(B.data() + base)); // 256-bits block from the base word.

    const auto unwanted_trail = i & 3;  // Number of unwanted bases (2-bits) trailing in the base word.
    if(!unwanted_trail)
        return blk;

    const auto to_clear_trail = _mm256_set1_epi64x(unwanted_trail * 2); // Number of trailing bits to clear from each word.
    const auto cleared = _mm256_srlv_epi64(blk, to_clear_trail);    // Trailing bits cleared from each 64-bit word.
    const auto r_shifted = _mm256_permute4x64_epi64(blk, 0b00'11'10'01);    // Words right-shifted by 1 word. The top word is don't-care.
    const auto to_clear_lead = _mm256_set1_epi64x((32 - unwanted_trail) * 2);    // Number of leading bits to clear from each word of the right-shifted block.
    const auto lost_bits = _mm256_sllv_epi64(r_shifted, to_clear_lead); // Bits lost due to the inability of whole register-wise right-shift during clearance of unwanted-bits.
    const auto restored = _mm256_or_si256(cleared, lost_bits);  // Restored lost trailing bits from words 1, 2, and 3.

    return restored;
}

}



#endif

