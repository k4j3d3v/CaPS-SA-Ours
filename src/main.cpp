
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

template <typename InputT>
int construct_and_dump_sa_helper(std::vector<InputT>& text, const std::string& op_path, const std::string& ext_mem_path, size_t subproblem_count, size_t max_context) {
    const bool ext_mem = true;  // TODO: take input.
    std::ofstream output(op_path);
    std::size_t n = text.size();
    std::cerr << "Text length: " << n << ".\n";
    if(n <= std::numeric_limits<uint32_t>::max())
    {
        CaPS_SA::Suffix_Array<InputT, uint32_t> suf_arr(text.data(), text.size(), ext_mem, ext_mem_path, subproblem_count, max_context);
        ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
        // suf_arr.dump(output);
    }
    else
    {
        CaPS_SA::Suffix_Array<InputT, uint64_t> suf_arr(text.data(), text.size(), ext_mem, ext_mem_path, subproblem_count, max_context);
        ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
        // suf_arr.dump(output);
    }

    output.close();
    return 0;
}

int construct_and_dump_sa(std::string input_t, const std::string& ip_path, const std::string& op_path, const std::string& ext_mem_path, size_t subproblem_count, size_t max_context) {
    if (input_t == "t"){
        std::vector<char> text;
        read_input<char>(ip_path, text);
        construct_and_dump_sa_helper<char>(text, op_path, ext_mem_path, subproblem_count, max_context);
    } else {
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
            construct_and_dump_sa_helper<uint64_t>(text, op_path, ext_mem_path, subproblem_count, max_context);
            input.close();
        } else {
            std::vector<uint32_t> text;
            text.resize(length);
            input.read(reinterpret_cast<char*>(text.data()), length * sizeof(uint32_t));
            construct_and_dump_sa_helper<uint32_t>(text, op_path, ext_mem_path, subproblem_count, max_context);
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


