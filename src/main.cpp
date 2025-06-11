
#include "Suffix_Array.hpp"
#include "Genomic_Text.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <type_traits>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <iostream>
#include <cassert>


template <typename T_seq_, typename T_idx_>
void pretty_print(const CaPS_SA::Suffix_Array<T_seq_, T_idx_>& suf_arr, std::ofstream& output)
{
    const std::size_t n = suf_arr.n();
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.SA()[i] << " \n"[i == n - 1];
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.LCP()[i] << " \n"[i == n - 1];
}

template <typename T_seq_>
int construct_and_dump_sa_helper(std::vector<T_seq_>& text, const std::string& op_path, const std::string& ext_mem_path, const size_t subproblem_count, const size_t max_context, const bool genomic)
{
    const bool ext_mem = true;  // TODO: take input.
    constexpr T_seq_ sentinel = std::is_same<T_seq_, char>::value ? '$' : std::numeric_limits<T_seq_>::max();

    // text.pop_back();
    std::size_t n = text.size();
    std::cerr << "Text length: " << n << ".\n";

    for(std::size_t i = 0; i < 7; ++i)
        text.push_back(sentinel);

    std::ofstream output(op_path);
    const auto construct = [&](auto sz)
    {
        if(!genomic)
        {
            CaPS_SA::Suffix_Array<T_seq_, decltype(sz)> suf_arr(text.data(), sz, ext_mem, ext_mem_path, subproblem_count, max_context);
            ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
            suf_arr.dump(output);
        }
        else
        {
            assert((std::is_same<T_seq_, char>::value));

            const CaPS_SA::Genomic_Text G(reinterpret_cast<const char*>(text.data()), sz);
            CaPS_SA::Suffix_Array<CaPS_SA::Genomic_Text, decltype(sz)> suf_arr(&G, sz, ext_mem, ext_mem_path, subproblem_count, max_context);
            ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
            suf_arr.dump(output);
        }
    };

    n <= std::numeric_limits<uint32_t>::max() ? construct(static_cast<uint32_t>(n)) : construct(static_cast<uint64_t>(n));

    output.close();
    return 0;
}

int construct_and_dump_sa(std::string input_t, const std::string& ip_path, const std::string& op_path, const std::string& ext_mem_path, size_t subproblem_count, size_t max_context)
{
    if(input_t == "t" || input_t == "g")
    {
        std::vector<char> text;
        CaPS_SA::read_input<char>(ip_path, text);
        construct_and_dump_sa_helper<char>(text, op_path, ext_mem_path, subproblem_count, max_context, input_t == "g");
    }
    else
    {
        std::ifstream input(ip_path);
        if (!input) {
            std::cerr << ip_path << " : could not be opened\n";
            std::exit(EXIT_FAILURE);
        }
        uint64_t length;
        uint64_t max_char;
        input.read(reinterpret_cast<char*>(&length), sizeof(length));
        input.read(reinterpret_cast<char*>(&max_char), sizeof(max_char));

        if (max_char >= std::numeric_limits<int32_t>::max()) {
            std::vector<uint64_t> text;
            text.resize(length);
            input.read(reinterpret_cast<char*>(text.data()), length * sizeof(uint64_t));
            construct_and_dump_sa_helper<uint64_t>(text, op_path, ext_mem_path, subproblem_count, max_context, false);
            input.close();
        } else {
            std::vector<uint32_t> text;
            text.resize(length);
            input.read(reinterpret_cast<char*>(text.data()), length * sizeof(uint32_t));
            construct_and_dump_sa_helper<uint32_t>(text, op_path, ext_mem_path, subproblem_count, max_context, false);
            input.close();
        }
    }
    return 0;
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
        std::cerr << "Usage: CaPS_SA <input_path> <output_path> <work_path_prefix> <(optional) input type [default: 't']> <(optional)-subproblem-count> <(optional)-bounded-context>>\n";
        std::exit(EXIT_FAILURE);
    }

    const std::string ip_path(argv[1]);
    const std::string op_path(argv[2]);
    const std::string ext_mem_path(argv[3]);
    const std::string data_type(argc >= 5 ? argv[4] : "t");
    const std::size_t subproblem_count(argc >= 6 ? std::atoi(argv[4]) : 0);
    const std::size_t max_context(argc >= 7 ? std::atoi(argv[5]) : 0);

    return construct_and_dump_sa(data_type, ip_path, op_path, ext_mem_path, subproblem_count, max_context);
}


