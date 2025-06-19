
#ifndef CAPS_SA_GENOMIC_TEXT_HPP
#define CAPS_SA_GENOMIC_TEXT_HPP


#ifdef USE_SIMDE
    #include "x86/avx2.h"
#else
    #include <immintrin.h>
#endif

#include <cstdint>
#include <cstddef>
#include <cstring>
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
#ifndef USE_SIMDE
    __m256i load(std::size_t i) const;
#else
    simde__m256i load(std::size_t i) const;
#endif

    // Returns the 29-nucleobase block (8 bytes) from onward the `i`'th
    // nucleobase, in 64-bits little-endian. The highest 6-bits are zeroed.
    uint64_t load_word(std::size_t i) const;

    // Returns the LCP length of the `124 x N`-nucleobases prefix of the
    // suffixes `x` and `y`.
    template <std::size_t N> std::size_t LCP_unrolled(std::size_t x, std::size_t y) const;

    // Returns the LCP length of the suffixes `x` and `y`, with context-length
    // `ctx`. `N x 124` nucleobases of prefix comparisons are loop-unrolled.
    template <std::size_t N> std::size_t LCP(std::size_t x, std::size_t y, std::size_t ctx) const;

public:

    // Constructs a 2-bit-packed representation of the genomic text `T` of
    // length `n`.
    Genomic_Text(const char* T, std::size_t n);

    Genomic_Text(const Genomic_Text&)  = delete;
    Genomic_Text(Genomic_Text&&) = delete;
    Genomic_Text& operator=(const Genomic_Text&) = delete;
    Genomic_Text& operator=(Genomic_Text&&) = delete;

    // Returns the code of the nucleobase at index `idx` of the original text.
    uint8_t operator[](std::size_t idx) const;

    // Returns the LCP-length of the suffixes at indices `x` and `y`, with
    // context-length `ctx`.
    std::size_t LCP(std::size_t x, std::size_t y, std::size_t ctx) const;

    // Returns the LCP length of the suffixes `x` and `y`, with context-length
    // `ctx`.
    std::size_t lcp_unvectorized(std::size_t x, std::size_t y, std::size_t ctx) const;
};


inline uint8_t Genomic_Text::operator[](const std::size_t idx) const
{
    assert(idx < n_);

    const auto byte_idx = idx >> 2;
    const auto bit_idx = (idx & 0b11) * 2;

    return  (B[byte_idx] & (0b11 << bit_idx)) >> bit_idx;
}


#ifndef USE_SIMDE
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
#else
    inline simde__m256i Genomic_Text::load(const std::size_t i) const
    {
        assert(i + 124 <= n_);

        const auto base = i / 4;    // Base word's index.
        const auto blk = simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i*>(B.data() + base)); // 256-bits block from the base word.

        const auto unwanted_trail = i & 3;  // Number of unwanted bases (2-bits) trailing in the base word.
        if(!unwanted_trail)
            return blk;

        const auto to_clear_trail = simde_mm256_set1_epi64x(unwanted_trail * 2); // Number of trailing bits to clear from each word.
        const auto cleared = simde_mm256_srlv_epi64(blk, to_clear_trail);    // Trailing bits cleared from each 64-bit word.
        const auto r_shifted = simde_mm256_permute4x64_epi64(blk, 0b00'11'10'01);    // Words right-shifted by 1 word. The top word is don't-care.
        const auto to_clear_lead = simde_mm256_set1_epi64x((32 - unwanted_trail) * 2);    // Number of leading bits to clear from each word of the right-shifted block.
        const auto lost_bits = simde_mm256_sllv_epi64(r_shifted, to_clear_lead); // Bits lost due to the inability of whole register-wise right-shift during clearance of unwanted-bits.
        const auto restored = simde_mm256_or_si256(cleared, lost_bits);  // Restored lost trailing bits from words 1, 2, and 3.

        return restored;
    }
#endif

inline uint64_t Genomic_Text::load_word(const std::size_t i) const
{
    assert(i + 29 <= n_);

    const auto base = i >> 2;   // Word's index.
    const auto unwanted_trail = i & 3;  // Number of unwanted bases (2-bits) trailing in the word.

    uint64_t w;
    std::memcpy(static_cast<void*>(&w), B.data() + base, 8);
    return (w >> (unwanted_trail * 2)) & 0x03FF'FFFF'FFFF'FFFFull;
}


#ifndef USE_SIMDE
    template <std::size_t N>
    inline std::size_t Genomic_Text::LCP_unrolled(const std::size_t x, const std::size_t y) const
    {
        if constexpr(N == 0)
            return 0;
        else
        {
            const auto X = load(x);
            const auto Y = load(y);

            const auto eq_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(X, Y));
            const auto neq_mask = ~eq_mask & 0x7FFF'FFFF;   // Top byte is degenerate in a block.
            if(neq_mask)
            {
                auto const X_b = reinterpret_cast<const unsigned char*>(&X);
                auto const Y_b = reinterpret_cast<const unsigned char*>(&Y);
                const auto bytes_eq = __builtin_ctz(neq_mask);
                assert(X_b[bytes_eq] != Y_b[bytes_eq]);

                const auto bits_eq = (bytes_eq << 3) + (__builtin_ctz(X_b[bytes_eq] ^ Y_b[bytes_eq]));
                return bits_eq >> 1;
            }

            return 124 + LCP_unrolled<N - 1>(x + 124, y + 124);
        }
    }
#else
    template <std::size_t N>
    inline std::size_t Genomic_Text::LCP_unrolled(const std::size_t x, const std::size_t y) const
    {
        if constexpr(N == 0)
            return 0;
        else
        {
            const auto X = load(x);
            const auto Y = load(y);

            const auto eq_mask = simde_mm256_movemask_epi8(simde_mm256_cmpeq_epi8(X, Y));
            const auto neq_mask = ~eq_mask & 0x7FFF'FFFF;   // Top byte is degenerate in a block.
            if(neq_mask)
            {
                auto const X_b = reinterpret_cast<const unsigned char*>(&X);
                auto const Y_b = reinterpret_cast<const unsigned char*>(&Y);
                const auto bytes_eq = __builtin_ctz(neq_mask);
                assert(X_b[bytes_eq] != Y_b[bytes_eq]);

                const auto bits_eq = (bytes_eq << 3) + (__builtin_ctz(X_b[bytes_eq] ^ Y_b[bytes_eq]));
                return bits_eq >> 1;
            }

            return 124 + LCP_unrolled<N - 1>(x + 124, y + 124);
        }
    }
#endif

template <std::size_t N>
inline std::size_t Genomic_Text::LCP(const std::size_t x, const std::size_t y, const std::size_t ctx) const
{
    std::size_t lcp = 0;

    if constexpr(N == 0)    // TODO: this branch needs to be optimized for bounded-ctx SAs.
    {
        for(; lcp < ctx; ++lcp)
            if((*this)[x + lcp] != (*this)[y + lcp])
                break;

        return lcp;
    }
    else
    {
        while((ctx - lcp) >= N * 124)
        {
            const auto l = LCP_unrolled<N>(x + lcp, y + lcp);
            lcp += l;
            if(l < N * 124)
                return lcp;
        }

        return lcp + LCP<0>(x + lcp, y + lcp, ctx - lcp);
    }
}


inline std::size_t Genomic_Text::LCP(const std::size_t x, const std::size_t y, const std::size_t ctx) const
{
    assert(x + ctx <= n_ && y + ctx <= n_);

    if(ctx < 29)
        return LCP<0>(x, y, ctx);

    const auto w_x = load_word(x);
    const auto w_y = load_word(y);
    assert((w_x >> 58) == 0 && (w_y >> 58) == 0);

    w_x != w_y ?
        assert((static_cast<std::size_t>(__builtin_ctzll(w_x ^ w_y)) >> 1) == LCP<1>(x, y, ctx)) :
        assert(LCP<1>(x, y, ctx) >= 29);

    return w_x != w_y ?
            static_cast<std::size_t>(__builtin_ctzll(w_x ^ w_y)) >> 1 :
            29 + LCP<1>(x + 29, y + 29, ctx - 29);
}


inline std::size_t Genomic_Text::lcp_unvectorized(const std::size_t x, const std::size_t y, const std::size_t ctx) const
{
    return LCP<0>(x, y, ctx);
}

}



#endif

