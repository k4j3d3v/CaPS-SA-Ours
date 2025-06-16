
#include "Suffix_Array.hpp"
#include "Genomic_Text.hpp"
#include "utility.hpp"
#include "CLI11.hpp"

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
int construct_and_dump_sa_helper(
        std::vector<T_seq_>& text, 
        const std::string& op_path, 
        const std::string& ext_mem_path, 
        const size_t subproblem_count, 
        const size_t max_context, 
        const bool genomic, 
        const bool ext_mem, 
        const bool output_lcp, 
        const bool collate_extmem_result)
{
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
            CaPS_SA::Suffix_Array<T_seq_, decltype(sz)> suf_arr(text.data(), sz, ext_mem, ext_mem_path, subproblem_count, max_context, output_lcp);
            ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
            if (ext_mem and collate_extmem_result) {
                suf_arr.dump(output);
                suf_arr.remove_extmem_partitions();
            } else {
                suf_arr.dump(output);
            }
        }
        else
        {
            assert((std::is_same<T_seq_, char>::value));

            const CaPS_SA::Genomic_Text G(reinterpret_cast<const char*>(text.data()), sz);
            CaPS_SA::Suffix_Array<CaPS_SA::Genomic_Text, decltype(sz)> suf_arr(&G, sz, ext_mem, ext_mem_path, subproblem_count, max_context, output_lcp);
            ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();
            if (ext_mem and collate_extmem_result) {
                suf_arr.dump(output);
                suf_arr.remove_extmem_partitions();
            } else {
                suf_arr.dump(output);
            }
        }
    };

    n <= std::numeric_limits<uint32_t>::max() ? construct(static_cast<uint32_t>(n)) : construct(static_cast<uint64_t>(n));

    output.close();
    return 0;
}

int construct_and_dump_sa(
        std::string input_t, 
        const std::string& ip_path, 
        const std::string& op_path, 
        const std::string& ext_mem_path, 
        size_t subproblem_count, 
        size_t max_context, 
        const bool ext_mem, 
        const bool output_lcp, 
        const bool collate_extmem_result)
{
    if(input_t == "t" || input_t == "g")
    {
        std::vector<char> text;
        CaPS_SA::read_input<char>(ip_path, text);
        construct_and_dump_sa_helper<char>(text, op_path, ext_mem_path, subproblem_count, max_context, input_t == "g", ext_mem, output_lcp, collate_extmem_result);
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
            construct_and_dump_sa_helper<uint64_t>(text, op_path, ext_mem_path, subproblem_count, max_context, false, ext_mem, output_lcp, collate_extmem_result);
            input.close();
        } else {
            std::vector<uint32_t> text;
            text.resize(length);
            input.read(reinterpret_cast<char*>(text.data()), length * sizeof(uint32_t));
            construct_and_dump_sa_helper<uint32_t>(text, op_path, ext_mem_path, subproblem_count, max_context, false, ext_mem, output_lcp, collate_extmem_result);
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
    CLI::App app{"CaPS-SA driver"};
    argv = app.ensure_utf8(argv);

    std::string ip_path= "";
    app.add_option("input", ip_path, "input path")->required();

    std::string op_path = "";
    app.add_option("output", op_path, "output path")->required();

    std::string data_type = "t";
    app.add_option("--data-type", data_type, "type of input data [text: \"t\", genomic: \"g\", or integer: \"i\"]")->check( [](const std::string &s) -> std::string {
        if (s == "t" or s == "g" or s == "i") {
            return "";
        } else {
            return std::string("The provided argument to --data-type is invalid, it must be one of t, g, or i.");
        }
    });

    bool ext_mem = false;
    auto ext_mem_flag = app.add_flag("--ext-mem", ext_mem, "pass this flag to use external memor construction");

    bool output_lcp = false;
    app.add_flag("--output-lcp", output_lcp, "pass this flag to output the LCP array along with the SA");
 
    bool collate_extmem_result = false;
    app.add_flag("--collate-extmem-result", collate_extmem_result, "collate the external memory buckets into a single file")->needs(ext_mem_flag);
    
    std::size_t subproblem_count = 0;
    app.add_option("--subproblem-count", subproblem_count, "subproblem count to use");

    std::size_t max_context = 0;
    app.add_option("--bounded-context", max_context, "bounded context to use (default: unlimited)");

    CLI11_PARSE(app, argc, argv);
   
    std::string ext_mem_prefix = ext_mem ? op_path : "";
    return construct_and_dump_sa(data_type, ip_path, op_path, ext_mem_prefix, subproblem_count, max_context, ext_mem, output_lcp, collate_extmem_result);
}


