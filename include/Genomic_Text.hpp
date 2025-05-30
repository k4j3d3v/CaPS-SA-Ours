
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

}



#endif

