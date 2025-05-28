
#ifndef CAPS_SA_GENOMIC_TEXT_HPP
#define CAPS_SA_GENOMIC_TEXT_HPP



#include <cstdint>
#include <cstddef>

#include <vector>


namespace CaPS_SA
{

// Class to to represent genomic texts in a 2-bit-packed manner.
class Genomic_Text
{
    const std::size_t n_;   // Length of the input text.

    const std::size_t pack_sz;  // Size of the bit-packed text.
    std::vector<uint8_t> B; // The bit-packed text.

public:

    // Constructs a 2-bit-packed representation of the genomic text `T` of
    // length `n`.
    Genomic_Text(const char* T, std::size_t n);
};

}



#endif

