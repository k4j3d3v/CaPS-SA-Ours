
#include "Genomic_Text.hpp"
#include "parlay/parallel.h"

#include <cassert>


namespace CaPS_SA
{

Genomic_Text::Genomic_Text(const char* const T, const std::size_t n):
      n_(n)
    , pack_sz((n_ + 3) / 4)
    , B(pack_sz + 1)    // +1 byte to facilitate one fell swoop loading of the last 31 bytes of the bit-packed text.
{
    const auto code = [](char ch)
    {
        ch &= static_cast<char>(~32);
        assert(ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T');

        return ((ch >> 2) ^ (ch >> 1)) & 0b11;
    };

    const auto pack = [&](const std::size_t i)
    {
        assert(i < pack_sz);
        assert(4 * i + 3 < n_);

        B[i] =    (code(T[i << 2]) << 0)
                | (code(T[(i << 2) | 0b01]) << 2)
                | (code(T[(i << 2) | 0b10]) << 4)
                | (code(T[(i << 2) | 0b11]) << 6);
    };

    parlay::parallel_for(0, n_ / 4, pack);
    if(n_ % 4)
    {
        B[n_ / 4] = 0;
        for(std::size_t i = 0; i < (n_ % 4); ++i)
            B[n_ / 4] |= (code(T[(n_ / 4) * 4 + i]) << (2 * i));
    }
}

}
