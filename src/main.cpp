
#include "Suffix_Array.hpp"

#include <vector>
#include <string>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <iostream>
#include <filesystem>


template <typename T_seq_>
void read_input(const std::string& ip_path, std::vector<T_seq_>& text)
{
    std::error_code ec;
    const auto file_size = std::filesystem::file_size(ip_path, ec);

    if(ec)
    {
        std::cerr << ip_path << " : " << ec.message() << "\n";
        std::exit(EXIT_FAILURE);
    }

    assert(file_size % sizeof(T_seq_) == 0);

    text.resize(file_size / sizeof(T_seq_));
    std::ifstream input(ip_path);
    input.read(reinterpret_cast<char*>(text.data()), file_size);
    input.close();
}


template <typename T_seq_, typename T_idx_>
void pretty_print(const CaPS_SA::Suffix_Array<T_seq_, T_idx_>& suf_arr, std::ofstream& output)
{
    const std::size_t n = suf_arr.n();
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.SA()[i] << " \n"[i == n - 1];
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.LCP()[i] << " \n"[i == n - 1];
}


int main(int argc, char* argv[])
{
#ifndef NDEBUG
    std::cout << "Warning: Executing in Debug Mode.\n";
#endif

    // TODO: standardize the API.
    constexpr auto arg_count = 4;
    if(argc < arg_count)
    {
        std::cerr << "Usage: CaPS_SA <input_path> <output_path> <work_path_prefix> <(optional)-subproblem-count> <(optional)-bounded-context>>\n";
        std::exit(EXIT_FAILURE);
    }


    const std::string ip_path(argv[1]);
    const std::string op_path(argv[2]);
    const std::string ext_mem_path(argv[3]);
    const std::size_t subproblem_count(argc >= 5 ? std::atoi(argv[4]) : 0);
    const std::size_t max_context(argc >= 6 ? std::atoi(argv[5]) : 0);
    const bool ext_mem = true;  // TODO: take input.

    typedef char T_seq_;

    std::vector<T_seq_> text;
    read_input(ip_path, text);
/*
    constexpr char lookup[4] = {'A', 'C', 'T', 'G'};
    size_t len = text.size();
    parlay::blocked_for(0, text.size(), 65536, 
      [&, len](size_t i, size_t start, size_t end) {
        (void)i;
        for (size_t j = start; j < std::min(end, len); ++j) {
          char c = text[j];
          text[j] = lookup[((std::toupper(c) & 0x6) >> 1)];
        };
    });
*/
    std::ofstream output(op_path);

    std::size_t n = text.size();
    std::cerr << "Text length: " << n << ".\n";
    if(n <= std::numeric_limits<uint32_t>::max())
    {
        CaPS_SA::Suffix_Array<T_seq_, uint32_t> suf_arr(text.data(), text.size(), ext_mem, ext_mem_path, subproblem_count, max_context);
        ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
        // suf_arr.dump(output);
    }
    else
    {
        CaPS_SA::Suffix_Array<T_seq_, uint64_t> suf_arr(text.data(), text.size(), ext_mem, ext_mem_path, subproblem_count, max_context);
        ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
        // suf_arr.dump(output);
    }

    output.close();


    return 0;
}
