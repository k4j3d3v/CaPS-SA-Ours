
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
        output << "SA[" << i << "]: " << suf_arr.SA()[i]<<"= "<< suf_arr.text()[suf_arr.SA()[i]] << " \n"[i == n - 1];
    // for(std::size_t i = 0; i < n; ++i)
    //     output << suf_arr.LCP()[i] << " \n"[i == n - 1];
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

    const auto construct = [&](const auto& sa)
    {
        typedef std::remove_reference_t<decltype(sa)> const_var;
        typedef std::remove_const_t<const_var> var;

        auto& suf_arr = *const_cast<var*>(&sa);
        ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();

        const std::string lcp_path = op_path + ".lcp";
        std::ofstream output(op_path);
        std::ofstream output_lcp(lcp_path);
        if (ext_mem and collate_extmem_result)
        {
            suf_arr.dump(output);
            suf_arr.remove_extmem_partitions();
        }
        else
            suf_arr.dump_separate(output, output_lcp);

        // pretty_print(suf_arr, std::cout);
        output.close();
        output_lcp.close();
    };

    using namespace CaPS_SA;
    if(!genomic)
    {
        std::cout<<"Constructing SA for integer text\n";
        if(n <= std::numeric_limits<uint32_t>::max())
        {
            std::cout<<"Using 32-bit SA\n";
            construct(Suffix_Array<T_seq_, uint32_t>(text.data(), n, ext_mem, ext_mem_path, subproblem_count, max_context, output_lcp));
        }
        else
        {
            std::cout<<"Using 64-bit SA\n"; 
            construct(Suffix_Array<T_seq_, uint64_t>(text.data(), n, ext_mem, ext_mem_path, subproblem_count, max_context, output_lcp));
        }
    }
    else
    {
        assert((std::is_same<T_seq_, char>::value));

        const Genomic_Text G(reinterpret_cast<const char*>(text.data()), n);
        n <= std::numeric_limits<uint32_t>::max() ?
            construct(Suffix_Array<Genomic_Text, uint32_t>(&G, n, ext_mem, ext_mem_path, subproblem_count, max_context, output_lcp)) :
            construct(Suffix_Array<Genomic_Text, uint64_t>(&G, n, ext_mem, ext_mem_path, subproblem_count, max_context, output_lcp));
    }

    return 0;
}
uint64_t stream_size(std::istream& is) {
    auto pos = is.tellg();
    is.seekg(0, std::ios::end);
    auto size = is.tellg();
    is.seekg(pos);
    return static_cast<uint64_t>(size);
}

int construct_and_dump_sa(
        std::string input_t, 
        std::string symbol_width, 
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
        // QUIIII
        std::vector<uint8_t> text;
        CaPS_SA::read_input<uint8_t>(ip_path, text);
        construct_and_dump_sa_helper<uint8_t>(text, op_path, ext_mem_path, subproblem_count, max_context, input_t == "g", ext_mem, output_lcp, collate_extmem_result);
    }
    else
    {
        std::ifstream input(ip_path);
        if (!input) {
            std::cerr << ip_path << " : could not be opened\n";
            std::exit(EXIT_FAILURE);
        }
    
        uint64_t length = stream_size(input);
        std::cerr << "Input file size: " << length << " bytes.\n";

        // uint64_t max_char;
        // input.read(reinterpret_cast<char*>(&length), sizeof(length));
        // input.read(reinterpret_cast<char*>(&max_char), sizeof(max_char));

        if (symbol_width == "64") {
            std::cerr << "Input data type: 64-bit integer.\n";
            std::vector<uint64_t> text;
            length /= sizeof(uint64_t);
            text.resize(length);
            input.read(reinterpret_cast<char*>(text.data()), length * sizeof(uint64_t));
            construct_and_dump_sa_helper<uint64_t>(text, op_path, ext_mem_path, subproblem_count, max_context, false, ext_mem, output_lcp, collate_extmem_result);
            input.close();
        } else {
            std::cerr << "Input data type: 32-bit integer.\n";

            std::vector<uint32_t> text;
            length /= sizeof(uint32_t);
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
#ifdef USE_AVX_512
    std::cerr << "Using AVX-512.\n";
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
    std::string symbol_width = "32";

    app.add_option(
        "--symbol-width",
        symbol_width,
        "Symbol width for integer inputs (32 or 64)"
    )->check(CLI::IsMember({"32", "64"}));

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

    std::size_t threads = 1;
    app.add_option("--threads", threads, "number of threads to use")->default_val("1");

    CLI11_PARSE(app, argc, argv);

    setenv("PARLAY_NUM_THREADS", std::to_string(threads).c_str(), 1);
   
    std::string ext_mem_prefix = ext_mem ? op_path : "";
    return construct_and_dump_sa(data_type, symbol_width, ip_path, op_path, ext_mem_prefix, subproblem_count, max_context, ext_mem, output_lcp, collate_extmem_result);
}


