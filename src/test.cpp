
#include "Suffix_Array.hpp"
#include "Genomic_Text.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <atomic>
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>
#include <fstream>
#include <filesystem>
#include <algorithm>


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
            const auto lcp_exp = lcp_unvectorized(T + i, T + j, n - std::max(i, j));
            const auto lcp_vec = suf_arr.LCP(i, j, n - std::max(i, j));
            const auto lcp_pck = G.LCP(i, j, n - std::max(i, j));

            if(lcp_vec != lcp_exp || lcp_pck != lcp_exp)
            {
                std::cerr << "At pair (" << i << ", " << j << "), expected LCP " << lcp_exp << "; LCP-vectorized " << lcp_vec << ", LCP-packed: " << lcp_pck << ".\n";
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
                bytes_scanned += (SA.LCP(i, j, n - j) + 1);

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
        std::cerr << "LCP-bandwidth on packed text: " << ((bases_scanned / duration(t_1 - t_0)) / (1024.0 * 1024.0)) << " Mbases/s.\n";
    }
// */
}


void LCP_bandwidth_thresholded(const char* const T, const std::size_t n, const std::size_t th = 32)
{
// /*
    {
        const Suffix_Array<char, uint32_t> SA(T, n);

        uint64_t bytes_scanned = 0;
        double time = 0;

        for(std::size_t i = 0; i < n; ++i)
            for(std::size_t j = i; j < n; ++j)
            {
                const auto t_pre = now();
                const auto scan = SA.LCP(i, j, n - j) + 1;
                if(scan >= th)
                {
                    const auto t_post = now();
                    bytes_scanned += scan;
                    time += duration(t_post - t_pre);
                }
            }

        std::cerr << "LCP-bandwidth on raw text with threshold " << th << ":    " << ((bytes_scanned / time) / (1024.0 * 1024.0)) << " MB/s.\n";
    }
// */

// /*
    {
        const Genomic_Text G(T, n);

        uint64_t bases_scanned = 0;
        double time = 0;

        for(std::size_t i = 0; i < n; ++i)
            for(std::size_t j = i; j < n; ++j)
            {
                const auto t_pre = now();
                const auto scan = G.LCP(i, j, n - j) + 1;
                if(scan >= th)
                {
                    const auto t_post = now();
                    bases_scanned += scan;
                    time += duration(t_post - t_pre);
                }
            }

        std::cerr << "LCP-bandwidth on packed text with threshold " << th << ": " << ((bases_scanned / time) / (1024.0 * 1024.0)) << " Mbases/s.\n";
    }
// */
}


template <typename T_idx_>
void LCP_bandwidth_selected(const char* const T, const std::size_t n, const std::string pairs_path)
{
    typedef std::pair<T_idx_, T_idx_> pair_t;

    constexpr uint64_t batch_sz = 256 * 1024 * 1024 / sizeof(std::pair<T_idx_, T_idx_>); // 256 MB.
    const uint64_t pairs_c = std::filesystem::file_size(pairs_path) / sizeof(pair_t);
    std::vector<pair_t> pairs(batch_sz);

    std::cerr << "#Pairs to compute LCP on: " << pairs_c << ".\n";

// /*
    {
        const Suffix_Array<char, T_idx_> SA(T, n);
        uint64_t bytes_scanned = 0;
        double time = 0;

        std::ifstream is(pairs_path, std::ios::binary);

        std::cerr << "\n";
        uint64_t pairs_read = 0;
        while(pairs_read < pairs_c)
        {
            const auto to_read = std::min(batch_sz, pairs_c - pairs_read);
            pairs.resize(to_read);
            is.read(reinterpret_cast<char*>(pairs.data()), to_read * sizeof(pair_t));
            pairs_read += to_read;

            const auto t_0 = CaPS_SA::now();

            for(const auto& p : pairs)
                bytes_scanned += (SA.LCP(p.first, p.second, n - p.second) + 1);

            const auto t_1 = CaPS_SA::now();
            time += duration(t_1 - t_0);

            std::cerr << "\rProcessed " << pairs_read << " pairs.";
        }

        is.close();

        std::cerr << "\n LCP-bandwidth on selected pairs on raw text:    " << ((bytes_scanned / time) / (1024.0 * 1024.0)) << " MB/s.\n";
    }
// */

// /*
    {
        const Genomic_Text G(T, n);
        uint64_t bases_scanned = 0;
        double time = 0;

        std::ifstream is(pairs_path, std::ios::binary);

        std::cerr << "\n";
        uint64_t pairs_read = 0;
        while(pairs_read < pairs_c)
        {
            const auto to_read = std::min(batch_sz, pairs_c - pairs_read);
            pairs.resize(to_read);
            is.read(reinterpret_cast<char*>(pairs.data()), to_read * sizeof(pair_t));
            pairs_read += to_read;

            const auto t_0 = CaPS_SA::now();

            for(const auto& p : pairs)
                bases_scanned += (G.LCP(p.first, p.second, n - p.second) + 1);

            const auto t_1 = CaPS_SA::now();
            time += duration(t_1 - t_0);

            std::cerr << "\rProcessed " << pairs_read << " pairs.";
        }

        is.close();

        std::cerr << "\n LCP-bandwidth on selected pairs on packed text: " << ((bases_scanned / time) / (1024.0 * 1024.0)) << " Mbases/s.\n";
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
    for(std::size_t i = 0; i < 7; ++i)
        text.push_back('$');

    if(n > std::numeric_limits<uint32_t>::max())
    {
        std::cerr << "Too large input for quadratic benchmark. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // CaPS_SA::cross_check_LCP(text.data(), n);
    // CaPS_SA::LCP_bandwidth(text.data(), n);
    // CaPS_SA::LCP_bandwidth_thresholded(text.data(), n, 32);

    if(argc >= 3)
    {
        const std::string pairs_path(argv[2]);
        CaPS_SA::LCP_bandwidth_selected<uint32_t>(text.data(), n, argv[2]);
    }

    return 0;
}
