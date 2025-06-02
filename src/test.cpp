
#include "Suffix_Array.hpp"
#include "Genomic_Text.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <atomic>
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <vector>


namespace CaPS_SA
{

void cross_check_LCP(const char* const T, const std::size_t n)
{
    const Suffix_Array<char, uint32_t> suf_arr(T, n);
    const Genomic_Text G(T, n);

    std::atomic_uint64_t solved = 0;
    parlay::parallel_for(0, n, [&](const std::size_t i)
    {
        parlay::parallel_for(0, n, [&](const std::size_t j)
        {
            const auto lcp_exp = suf_arr.LCP(T + i, T + j, n - std::max(i, j));
            const auto lcp_comp = G.LCP(i, j, n - std::max(i, j));

            if(lcp_comp != lcp_exp)
            {
                std::cerr << "At pair (" << i << ", " << j << "), expected LCP " << lcp_exp << ", received LCP " << lcp_comp << ".\n";
                std::exit(EXIT_FAILURE);
            }
        }
        , 1);

        solved++;

        const uint64_t s = solved;
        if(s % 1024 == 0)
            std::cerr << "\rChecked " << s << " starting positions.";
    });

    std::cerr << "\rChecked " << solved << " starting positions.\n";
}


void LCP_bandwidth(const char* const T, const std::size_t n)
{
// /*
    {
        const Suffix_Array<char, uint32_t> SA(T, n);

        const auto t_0 = CaPS_SA::now();

        uint64_t bytes_scanned = 0;
        for(std::size_t i = 0; i < n; ++i)
            for(std::size_t j = i; j < n; ++j)
                bytes_scanned += (SA.LCP(T + i, T + j, n - j) + 1);

        const auto t_1 = CaPS_SA::now();
        std::cerr << "LCP-bandwidth on raw text:    " << ((bytes_scanned / duration(t_1 - t_0)) / (1024.0 * 1024.0)) << " MB/s.\n";
    }
// */

// /*
    {
        const Genomic_Text G(T, n);

        const auto t_0 = CaPS_SA::now();

        uint64_t bases_scanned = 0;
        for(std::size_t i = 0; i < n; ++i)
            for(std::size_t j = i; j < n; ++j)
                bases_scanned += (G.LCP(i, j, n - j) + 1);

        const auto t_1 = CaPS_SA::now();
        std::cerr << "LCP-bandwidth on packed text: " << (((bases_scanned / 4.0) / duration(t_1 - t_0)) / (1024.0 * 1024.0)) << " MB/s.\n";
    }
// */
}

}


int main(int argc, char* argv[])
{
#ifndef NDEBUG
    std::cout << "Warning: Executing in Debug Mode.\n";
#endif

    (void)argc, (void)argv;

    const std::string ip_path(argv[1]);

    std::vector<char> text;
    CaPS_SA::read_input(ip_path, text);

    text.pop_back();
    for(char& ch : text)
        if(ch > 'T')
            ch = std::toupper(ch);

    const std::size_t n = text.size();

    if(n > std::numeric_limits<uint32_t>::max())
    {
        std::cerr << "Too large input for quadratic benchmark. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    CaPS_SA::cross_check_LCP(text.data(), n);
    CaPS_SA::LCP_bandwidth(text.data(), n);

    return 0;
}
