
#include "Suffix_Array.hpp"
#include "parlay/parallel.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace CaPS_SA
{

template <typename T_idx_>
Suffix_Array<T_idx_>::Suffix_Array(const char* const T, const idx_t n, const bool ext_mem, const std::string& ext_mem_path, const idx_t subproblem_count, const idx_t max_context):
    T_(T),
    n_(n),
    p_(subproblem_count > 0 ? subproblem_count : default_subproblem_count),
    per_worker_in_mem_elem(!ext_mem ? 0 : 2 * static_cast<idx_t>(std::ceil(n_ / p_))),
    SA_(!ext_mem ? allocate<idx_t>(n_) : nullptr),
    LCP_(!ext_mem ? allocate<idx_t>(n_) : nullptr),
    SA_w(nullptr),
    LCP_w(nullptr),
    ext_mem_(ext_mem),
    ext_mem_path(ext_mem_path),
    max_context(max_context ? max_context : n_),
    sample_per_part_(static_cast<idx_t>(std::ceil(32.0 * std::log(n_)))),   // c \ln n
    pivot_(nullptr),
    part_size_scan_(nullptr),
    part_ruler_(nullptr)
{
    if(p_ == 0)
    {
        std::cerr << "The environment variable `PARLAY_NUM_THREADS` needs to be set. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_idx_>
Suffix_Array<T_idx_>::~Suffix_Array()
{
    if(!ext_mem_)
        std::free(SA_),
        std::free(LCP_);
    else
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            std::free(SA_buf[w_id].data),
            std::free(LCP_buf[w_id].data),
            std::free(SA_w_buf[w_id].data),
            std::free(LCP_w_buf[w_id].data);
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const
{
    idx_t m = 0;    // LCP of the last compared pair.
    idx_t l_x;  // LCP(X_i, X_{i - 1}).
    idx_t i = 0;    // Index into `X`.
    idx_t j = 0;    // Index into `Y`.
    idx_t k = 0;    // Index into `Z`.

    while(i < len_x && j < len_y)
    {
        l_x = LCP_x[i];

        if(l_x > m)
            Z[k] = X[i],
            LCP_z[k] = l_x,
            m = m;
        else if(l_x < m)
            Z[k] = Y[j],
            LCP_z[k] = m,
            m = l_x;
        else    // Compute LCP of X_i and Y_j through linear scan.
        {
            const idx_t max_n = n_ - std::max(X[i], Y[j]);  // Length of the shorter suffix.
            const idx_t context = std::min(max_context, max_n); // Prefix-context length for the suffixes.
            const idx_t n = m + lcp_opt_avx(T_ + (X[i] + m), T_ + (Y[j] + m), context - m); // LCP(X_i, Y_j)

            // Whether the shorter suffix is a prefix of the longer one.
            Z[k] = (n == max_n ?    std::max(X[i], Y[j]) :
                                    (T_[X[i] + n] < T_[Y[j] + n] ? X[i] : Y[j]));
            LCP_z[k] = (Z[k] == X[i] ? l_x : m);
            m = n;
        }


        if(Z[k] == X[i])
            i++;
        else    // Swap X and Y (and associated data structures) when Y_j gets pulled into Z.
        {
            j++;
            std::swap(X, Y),
            std::swap(len_x, len_y),
            std::swap(LCP_x, LCP_y),
            std::swap(i, j);
        }

        k++;
    }


    for(; i < len_x; ++i, ++k)  // Copy rest of the data from X to Z.
        Z[k] = X[i], LCP_z[k] = LCP_x[i];

    if(j < len_y)   // Copy rest of the data from Y to Z.
    {
        Z[k] = Y[j], LCP_z[k] = m;
        for(j++, k++; j < len_y; ++j, ++k)
            Z[k] = Y[j], LCP_z[k] = LCP_y[j];
    }
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::merge_sort(idx_t* const X, idx_t* const Y, const idx_t n, idx_t* const LCP, idx_t* const W) const
{
    assert(std::memcmp(X, Y, n * sizeof(idx_t)) == 0);

    if(n == 1)
        LCP[0] = 0;
    else
    {
        const idx_t m = n / 2;
        const auto f = [&](){ merge_sort(Y, X, m, W, LCP); };
        const auto g = [&](){ merge_sort(Y + m, X + m, n - m, W + m, LCP + m); };

        m < nested_par_grain_size ?
            (f(), g()) : parlay::par_do(f, g, ext_mem_);
        merge(X, m, X + m, n - m, W, W + m, Y, LCP);
    }
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::initialize()
{
    const auto t_s = now();

    if(!ext_mem_)
    {
        SA_w = allocate<idx_t>(n_);
        LCP_w = allocate<idx_t>(n_);
    }
    else
    {
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            SA_buf.emplace_back(allocate<idx_t>(per_worker_in_mem_elem)),
            LCP_buf.emplace_back(allocate<idx_t>(per_worker_in_mem_elem)),
            SA_w_buf.emplace_back(allocate<idx_t>(per_worker_in_mem_elem)),
            LCP_w_buf.emplace_back(allocate<idx_t>(per_worker_in_mem_elem)),
            pivot_loc_buf.emplace_back(allocate<idx_t>(p_ + 2));

        for(std::size_t p_id = 0; p_id < p_; ++p_id)
            SA_bucket.emplace_back(ext_mem_path + "_SA_" + std::to_string(p_id)),
            LCP_bucket.emplace_back(ext_mem_path + "_LCP_" + std::to_string(p_id)),
            sz_bucket.emplace_back(ext_mem_path + "_sz_" + std::to_string(p_id));

        lock = new Padded_Data<Spin_Lock>[p_];
    }

    const auto sample_count = p_ * sample_per_part_;
    pivot_ = allocate<idx_t>(sample_count);

    const auto t_e = now();
    std::cerr << "Initialized required data structures. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sort_subarrays()
{
    const auto t_s = now();

    const auto mem_init = [SA_ = SA_, SA_w = SA_w](const idx_t i){ SA_[i] = SA_w[i] = i; };
    parlay::parallel_for(0, n_, mem_init);

    const auto subarr_size = n_ / p_;   // Size of each subarray to be sorted independently.
    const auto sort_subarr =
        [&](const idx_t i)
        {
            merge_sort( SA_w + i * subarr_size, SA_ + i * subarr_size,
                        subarr_size + (i < p_ - 1 ? 0 : n_ % p_),
                        LCP_ + i * subarr_size, LCP_w + i * subarr_size);

            if(++solved_ % 8 == 0)
                std::cerr << "\rSorted " << solved_ << " subarrays.";
        };

    solved_ = 0;
    parlay::parallel_for(0, p_, sort_subarr, 1);
    std::cerr << "\n";

    const auto t_e = now();
    std::cerr << "Sorted the subarrays independently. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sort_subarrays_ext_mem()
{
    const auto t_s = now();

    const auto subarr_sz = n_ / p_; // Size of each subarray to sort independently.
    const auto sort_subarr =
        [&](const idx_t p_id)
        {
            const auto w_id = parlay::worker_id();
            const auto range_beg = p_id * subarr_sz, len = subarr_sz + (p_id < p_ - 1 ? 0 : n_ % p_);
            assert(len < per_worker_in_mem_elem);

            const auto SA = SA_buf[w_id].data, SA_w = SA_w_buf[w_id].data, LCP = LCP_buf[w_id].data, LCP_w = LCP_w_buf[w_id].data;
            for(std::size_t i = 0; i < len; ++i)
                SA[i] = SA_w[i] = range_beg + i;

            merge_sort(SA_w, SA, len, LCP, LCP_w);
            assert(is_sorted(SA, len, LCP));

            if(++solved_ % 8 == 0)
                std::cerr << "\rSorted and partitioned " << solved_ << " subarrays.";

            // const auto pivot_off = p_id * sample_per_part_;
            // sample_pivots(SA, len, sample_per_part_, pivot_ + pivot_off);


            auto const P = pivot_loc_buf[w_id].data;    // Pivot positions in this sorted subarray.

            P[0] = 0, P[p_] = len;  // Two flanking pivot indices.
            for(idx_t piv_id = 0; piv_id < p_ - 1; ++piv_id)
                P[piv_id + 1] = upper_bound(SA, len, T_ + pivot_[piv_id], n_ - pivot_[piv_id]);

            distribute_sub_subarrays_ext_mem();
        };

    solved_ = 0;
    parlay::parallel_for(0, p_, sort_subarr, 1);
    std::cerr << "\n";

    for(idx_t p_id = 0; p_id < p_; ++p_id)
        SA_bucket[p_id].data.close(),
        LCP_bucket[p_id].data.close(),
        sz_bucket[p_id].data.close();

    const auto t_e = now();
    std::cerr << "Sorted the subarrays independently. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sample_pivots(const idx_t* const X, const idx_t n, const idx_t m, idx_t* const P)
{
    const auto gap = n / (m + 1);   // Distance-gap between pivots.
    for(idx_t i = 0; i < m; ++i)
        P[i] = X[(i + 1) * gap - 1];
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sample_pivots(const idx_t r_beg, const idx_t n, const idx_t m, idx_t* const P)
{
    const auto gap = n / (m + 1);
    for(idx_t i = 0; i < m; ++i)
        P[i] = r_beg + (i + 1) * gap - 1;
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::select_pivots()
{
    const auto t_s = now();

    !ext_mem_ ? collect_samples() : collect_samples_ext_mem();
    select_pivots_off_samples();

    const auto t_e = now();
    std::cerr << "Selected the global pivots. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::collect_samples()
{
    const auto subarr_size = n_ / p_;   // Size of each sorted subarray.
    for(idx_t i = 0; i < p_; ++i)
        sample_pivots(  SA_ + i * subarr_size, subarr_size + (i < p_ - 1 ? 0 : n_ % p_),
                        sample_per_part_, pivot_ + i * sample_per_part_);
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::collect_samples_ext_mem()
{
    const auto subarr_sz = n_ / p_; // Size of each sorted subarray.
    for(idx_t i = 0; i < p_; ++i)
        sample_pivots(  i * subarr_sz, subarr_sz + (i < p_ - 1? 0 : n_ % p_),
                        sample_per_part_, pivot_ + i * sample_per_part_);
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::select_pivots_off_samples()
{
    const auto sample_count = p_ * sample_per_part_;    // Total number of samples to select pivots from.
    idx_t* const pivot_w = allocate<idx_t>(sample_count);   // Working space to sample pivots.
    auto const temp_1 = allocate<idx_t>(sample_count), temp_2 = allocate<idx_t>(sample_count);

    std::memcpy(pivot_w, pivot_, sample_count * sizeof(idx_t));
    merge_sort(pivot_, pivot_w, sample_count, temp_1, temp_2);

    sample_pivots(pivot_w, sample_count, p_ - 1, pivot_);

    std::free(pivot_w), std::free(temp_1), std::free(temp_2);
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::locate_pivots(idx_t* const P) const
{
    const auto t_s = now();

    const auto subarr_size = n_ / p_;   // Size of each independent sorted subarray.

    const auto locate =
        [&](const idx_t i)
        {
            const auto X_i = SA_ + i * subarr_size; // The i'th subarray.
            const auto P_i = P + i * (p_ + 1);  // Pivot locations in `X_i` are to be placed in `P_i`.

            P_i[0] = 0, P_i[p_] = subarr_size + (i < p_ - 1 ? 0 : n_ % p_); // The two flanking pivot indices.
            for(idx_t j = 0; j < p_ - 1; ++j) // TODO: try parallelizing this loop too; observe performance diff.
                P_i[j + 1] = upper_bound(X_i, P_i[p_], T_ + pivot_[j], n_ - pivot_[j]);
        };

    parlay::parallel_for(0, p_, locate, 1);

    const auto t_e = now();
    std::cerr << "Located the pivots in each sorted subarray. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
T_idx_ Suffix_Array<T_idx_>::upper_bound(const idx_t* const X, const idx_t n, const char* const P, const idx_t P_len) const
{
    // Invariant: SA[l] < s < SA[r].

    int64_t l = -1, r = n;  // (Exclusive-) Range of the iterations in the binary search.
    idx_t c;    // Midpoint in each iteration.
    idx_t soln = n; // Solution of the search.
    idx_t lcp_l = 0, lcp_r = 0; // LCP(s, SA[l]) and LCP(s, SA[r]).
	idx_t approx = 65536;   // TODO: better tune and document.

    while(r - l > 1)    // Candidate matches exist.
    {
        c = (l + r) / 2;
        const char* const suf = T_ + X[c];  // The suffix at the middle.
        const auto suf_len = n_ - X[c]; // Length of the suffix.

        idx_t lcp_c = std::min(lcp_l, lcp_r);   // LCP(X[c], P).
        lcp_c = std::min(lcp_c, approx);   // LCP(X[c], P).
        auto max_lcp = std::min(suf_len, P_len);  // Maximum possible LCP, i.e. length of the shorter string.
		max_lcp = std::min(max_lcp,approx);
        lcp_c += lcp_opt_avx(suf + lcp_c, P + lcp_c, max_lcp - lcp_c);  // Skip an informed number of character comparisons.

        if(lcp_c == max_lcp)    // One is a prefix of the other.
        {
            if(lcp_c == P_len)  // P is a prefix of the suffix.
            {
                if(P_len == suf_len)    // The query is the suffix itself, i.e. P = X[c]
                    return c + 1;
                else    // P < X[c]
                    r = c, lcp_r = lcp_c, soln = c;
            }
            else    // The suffix is a prefix of the query, so X[c] < P; technically impossible if the text terminates with $.
                l = c, lcp_l = lcp_c;
        }
        else    // Neither is a prefix of the other.
            if(suf[lcp_c] < P[lcp_c])   // X[c] < P
                l = c, lcp_l = lcp_c;
            else    // P < X[c]
                r = c, lcp_r = lcp_c, soln = c;
    }


    return soln;
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::partition_sub_subarrays(const idx_t* const P)
{
    const auto t_s = now();

    part_size_scan_ = allocate<idx_t>(p_ + 1);

    const auto collect_size =   // Collects the size of the `j`'th partition.
        [&](const idx_t j)
        {
            part_size_scan_[j] = 0;
            for(idx_t i = 0; i < p_; ++i)   // For subarray `i`.
            {
                const auto P_i = P + i * (p_ + 1);  // Pivot collection of subarray `i`.
                part_size_scan_[j] += (P_i[j + 1] - P_i[j]);    // Collect its `j`'th sub-subarray's size.
            }
        };

    parlay::parallel_for(0, p_, collect_size, 1);   // Collect the individual size of each partition.


    // Compute inclusive-scan (prefix sum) of the partition sizes.
    idx_t curr_sum = 0;
    for(idx_t j = 0; j < p_; ++j) // For partition `j`.
    {
        const auto part_size = part_size_scan_[j];

        part_size_scan_[j] = curr_sum;
        curr_sum += part_size;
    }

    part_size_scan_[p_] = curr_sum;
    assert(part_size_scan_[p_] == n_);


    // Collate the sorted sub-subarrays to appropriate partitions.
    part_ruler_ = allocate<idx_t>(p_ * (p_ + 1));
    const idx_t subarr_size = n_ / p_;
    const auto collate =    // Collates the `j`'th sub-subarray from each sorted subarray to partition `j`.
        [&](const idx_t j)
        {
            auto const Y_j = SA_w + part_size_scan_[j]; // Memory-base for partition `j`.
            auto const LCP_Y_j = LCP_w + part_size_scan_[j];    // Memory-base for LCPs of partition `j`.
            auto const sub_subarr_idx = part_ruler_ + j * (p_ + 1); // Index of the sorted sub-subarrays in `Y_j`.
            idx_t curr_idx = 0; // Current index into `Y_j`.

            for(idx_t i = 0; i < p_; ++i)   // Subarray `i`.
            {
                const auto X_i = SA_ + i * subarr_size; // `i`'th sorted subarray.
                const auto LCP_X_i = LCP_ + i * subarr_size;    // LCP array of `X_i`.
                const auto P_i = P + i * (p_ + 1);  // Pivot collection of subarray `i`.

                const auto sub_subarr_size = P_i[j + 1] - P_i[j];   // Size of the `j`'th sub-subarray of subarray `i`.
                sub_subarr_idx[i] = curr_idx;
                std::memcpy(Y_j + sub_subarr_idx[i], X_i + P_i[j], sub_subarr_size * sizeof(idx_t));
                std::memcpy(LCP_Y_j + sub_subarr_idx[i], LCP_X_i + P_i[j], sub_subarr_size * sizeof(idx_t));
                LCP_Y_j[sub_subarr_idx[i]] = 0;
                curr_idx += sub_subarr_size;
            }

            sub_subarr_idx[p_] = curr_idx;
            assert(curr_idx == part_size_scan_[j + 1] - part_size_scan_[j]);
        };

    parlay::parallel_for(0, p_, collate, 1);

    const auto t_e = now();
    std::cerr << "Collated the sorted sub-subarrays into partitions. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::distribute_sub_subarrays_ext_mem()
{
    const auto w_id = parlay::worker_id();
    const auto SA = SA_buf[w_id].data, LCP = LCP_buf[w_id].data, P = pivot_loc_buf[w_id].data;


    // Different sequences of parts-copying (dispersion) by different workers to minimize lock-contention.

    const auto step = parlay::num_workers();
    for(std::size_t round = 0; round < step; ++round)
        for(auto part_id = (w_id + round) % step; part_id < p_; part_id += step)
        {
            const auto sub_subarr_sz = P[part_id + 1] - P[part_id];

            lock[part_id].data.lock();

            auto &SA_b = SA_bucket[part_id].data, &LCP_b = LCP_bucket[part_id].data, &sz_b = sz_bucket[part_id].data;
            SA_b.add(SA + P[part_id], sub_subarr_sz);
            LCP_b.add(LCP + P[part_id], sub_subarr_sz);
            sz_b.add(sub_subarr_sz);

            lock[part_id].data.unlock();
        }
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::merge_sub_subarrays()
{
    const auto t_s = now();

    const auto mem_init =
        [&](const idx_t j)
        {
            const auto part_size = part_size_scan_[j + 1] - part_size_scan_[j];
            std::memcpy(SA_ + part_size_scan_[j], SA_w + part_size_scan_[j], part_size * sizeof(idx_t));
            std::memcpy(LCP_ + part_size_scan_[j], LCP_w + part_size_scan_[j], part_size * sizeof(idx_t));
        };

    parlay::parallel_for(0, p_, mem_init, 1);   // Fulfill `sort_partition`'s precondition.


    const auto sort_part =
        [&](const idx_t j)
        {
            const auto part_idx = part_size_scan_[j];   // Index of the partition in the partitions' flat collection.
            auto const X_j = SA_w + part_idx;   // Memory-base for partition `j`.
            auto const Y_j = SA_ + part_idx;    // Location to sort partition `j`.
            auto const LCP_X_j = LCP_w + part_idx;  // Memory-base for the LCP-arrays of partition `j`.
            auto const LCP_Y_j = LCP_ + part_idx;   // LCP array of `Y_j`.
            auto const sub_subarr_idx = part_ruler_ + j * (p_ + 1); // Indices of the sorted subarrays in `X_i`.

            sort_partition(X_j, Y_j, p_, sub_subarr_idx, LCP_X_j, LCP_Y_j);

            if(++solved_ % 8 == 0)
                std::cerr << "\rMerged " << solved_ << " partitions.";
        };

    solved_ = 0;
    parlay::parallel_for(0, p_, sort_part, 1);  // Merge the sorted subarrays in each partitions.
    std::cerr << "\n";

    const auto t_e = now();
    std::cerr << "Merged the sorted subarrays in each partition. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sort_partition(idx_t* const X, idx_t* const Y, const idx_t n, const idx_t* const S, idx_t* const LCP_x, idx_t* const LCP_y)
{
    if(n == 1)
        return;

    const auto m = n / 2;
    const auto flat_count_l = S[m] - S[0];
    const auto flat_count_r = S[n] - S[m];

    const auto f = [&](){ sort_partition(Y, X, m, S, LCP_y, LCP_x); };
    const auto g = [&](){ sort_partition(Y + flat_count_l, X + flat_count_l, n - m, S + m, LCP_y + flat_count_l, LCP_x + flat_count_l); };

    (flat_count_l < nested_par_grain_size || flat_count_r < nested_par_grain_size) ?
        (f(), g()) : parlay::par_do(f, g);
    merge(X, flat_count_l, X + flat_count_l, flat_count_r, LCP_x, LCP_x + flat_count_l, Y, LCP_y);
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::compute_partition_boundary_lcp()
{
    const auto t_s = now();

    const auto compute_boundary_lcp =
        [&](const idx_t j)
        {
          const auto part_idx = part_size_scan_[j];
          LCP_[part_idx] = lcp_opt_avx(T_ + SA_[part_idx - 1], T_ + SA_[part_idx], n_ - std::max(SA_[part_idx - 1], SA_[part_idx]));
        };

    parlay::parallel_for(1, p_, compute_boundary_lcp, 1);

    const auto t_e = now();
    std::cerr << "Computed the LCPs at the partition boundaries. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::clean_up()
{
    const auto t_s = now();

    std::free(SA_w);
    std::free(LCP_w);

    std::free(pivot_);
    std::free(part_size_scan_);
    std::free(part_ruler_);

    const auto t_e = now();
    std::cerr << "Released the temporary data structures. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::construct()
{
    const auto t_start = now();

    initialize();

    // merge_sort(SA_w, SA_, n_, LCP_, LCP_w);  // Monolithic construction.
    sort_subarrays();

    select_pivots();

    idx_t* const P = allocate<idx_t>(p_ * (p_ + 1));  // Collection of pivot locations in the subarrays.
    locate_pivots(P);
    partition_sub_subarrays(P);
    std::free(P);

    merge_sub_subarrays();

    compute_partition_boundary_lcp();

    clean_up();

    const auto t_end = now();
    std::cerr << "Constructed the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::construct_ext_mem()
{
    const auto t_start = now();

    initialize();

    select_pivots();

    sort_subarrays_ext_mem();
    // select_pivots_off_samples();

    const auto t_end = now();
    std::cerr << "Constructed the suffix array and the LCP array. Time taken: " << duration(t_end - t_start) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::dump(std::ofstream& output)
{
    const auto t_start = now();

    const std::size_t n = n_;
    output.write(reinterpret_cast<const char*>(&n), sizeof(std::size_t));
    output.write(reinterpret_cast<const char*>(SA_), n_ * sizeof(idx_t));
    output.write(reinterpret_cast<const char*>(LCP_), n_ * sizeof(idx_t));

    const auto t_end = now();
    std::cerr << "Dumped the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n";
}


template <typename T_idx_>
bool Suffix_Array<T_idx_>::is_sorted(const idx_t* const X, const idx_t n, const idx_t* const LCP_X) const
{
    if(LCP_X && n > 0 && LCP_X[0] != 0)
        return false;

    for(idx_t i = 1; i < n; ++i)
    {
        const auto x = T_ + X[i - 1], y = T_ + X[i];
        const auto l = std::min(n_ - X[i - 1], n_ - X[i]);

        idx_t lcp = 0;
        for(; lcp < l; ++lcp)
            if(x[lcp] < y[lcp])
                break;
            else if(x[lcp] > y[lcp])
                return false;

        if(LCP_X && lcp != LCP_X[i])
            return false;
    }

    return true;
}

}



// Template instantiations for the required instances.
template class CaPS_SA::Suffix_Array<uint32_t>;
template class CaPS_SA::Suffix_Array<uint64_t>;
