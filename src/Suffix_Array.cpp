
#include "Suffix_Array.hpp"
#include "parlay/parallel.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <numeric>
#include <random>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace CaPS_SA
{

template <typename T_seq_, typename T_idx_>
Suffix_Array<T_seq_, T_idx_>::Suffix_Array(const T_seq_* const T, const idx_t n, const bool ext_mem, const std::string& ext_mem_path, const idx_t subproblem_count, const idx_t max_context, const bool output_lcp):
    T_(T),
    n_(n),
    p_(std::min(subproblem_count > 0 ? subproblem_count : default_subproblem_count, n / 16)),   // TODO: fix subproblem-count for small `n`.
    per_worker_in_mem_elem(!ext_mem ? 0 : static_cast<idx_t>(std::ceil(n_ / p_))),
    SA_(!ext_mem ? allocate<idx_t>(n_) : nullptr),
    LCP_(!ext_mem ? allocate<idx_t>(n_) : nullptr),
    SA_w(nullptr),
    LCP_w(nullptr),
    ext_mem_ctr_(ext_mem),
    ext_mem_path(ext_mem_path),
    max_context(max_context ? max_context : n_),
    sample_per_part_(std::min(static_cast<idx_t>(std::ceil(32.0 * std::log(n_))), n_ / p_ - 1)),    // (c \ln n) or (|subarray| - 1)
    pivot_(nullptr),
    part_size_scan_(nullptr),
    part_ruler_(nullptr),
    op_lcp(output_lcp),   
    removed_lcp_partitions(false),
    removed_extmem_partitions(false),
    constructed(false)
{
    assert(n_ >= 16);   // TODO: fix subproblem-count for small `n`.

    if(p_ > n_)
    {
        CAPS_SA_LOG(std::cerr << "Incompatible subproblem-count. Aborting.\n");
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_seq_, typename T_idx_>
Suffix_Array<T_seq_, T_idx_>::~Suffix_Array()
{
    if(!ext_mem_ctr_)
        deallocate(SA_),
        deallocate(LCP_);
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const
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

            assert((X[i] + m) + (context - m) <= n_ && (Y[j] + m) + (context - m) <= n_);

            const idx_t n = m + LCP(X[i] + m, Y[j] + m, context - m);   // LCP(X_i, Y_j)

            // Whether the shorter suffix is a prefix of the longer one.
            Z[k] = (n == max_n ?    std::max(X[i], Y[j]) :
                                    (is_lesser(X[i] + n, Y[j] + n) ? X[i] : Y[j]));
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


    // Copy rest of the data from X to Z.
    std::memcpy(Z + k, X + i, (len_x - i) * sizeof(idx_t));
    std::memcpy(LCP_z + k, LCP_x + i, (len_x - i) * sizeof(idx_t));

    // Copy rest of the data from Y to Z.
    std::memcpy(Z + k, Y + j, (len_y - j) * sizeof(idx_t));
    std::memcpy(LCP_z + k, LCP_y + j, (len_y - j) * sizeof(idx_t));

    // Fix first copied LCP.
    if(k < len_x + len_y)   // Both collections should not be empty.
        LCP_z[k] = m;
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::merge_sort(idx_t* const X, idx_t* const Y, const idx_t n, idx_t* const LCP, idx_t* const W) const
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
            (f(), g()) : parlay::par_do(f, g, ext_mem_ctr_);
        merge(X, m, X + m, n - m, W, W + m, Y, LCP);
    }
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::initialize()
{
    const auto t_s = now();

    if(!ext_mem_ctr_)
    {
        SA_w = allocate<idx_t>(n_);
        LCP_w = allocate<idx_t>(n_);
        part_size_scan_ = allocate<idx_t>(p_ + 1);
        part_ruler_ = allocate<idx_t>(p_ * (p_ + 1));
    }
    else
    {
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            w_buf.emplace_back(*this);

        for(std::size_t p_id = 0; p_id < p_; ++p_id)
            subproblem_space.emplace_back(Subproblem_Ext_Mem(*this, p_id));

        std::vector<Padded<Spin_Lock>> temp(p_);
        lock.swap(temp);
    }

    const auto sample_count = p_ * sample_per_part_;
    pivot_ = allocate<idx_t>(sample_count);

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Initialized required data structures. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::permute()
{
    const auto t_s = now();

    const auto populate = [&](const idx_t i){ SA_[i] = SA_w[i] = i; };
    parlay::parallel_for(0, n_, populate);

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Populated the suffix array with a permutation. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::sort_subarrays()
{
    const auto t_s = now();

    const auto subarr_size = n_ / p_;   // Size of each subarray to be sorted independently.
    std::atomic_uint64_t solved_ = 0;   // Progress tracker—number of subproblems solved in some step.
    const auto sort_subarr =
        [&](const idx_t i)
        {
            merge_sort( SA_w + i * subarr_size, SA_ + i * subarr_size,
                        subarr_size + (i < p_ - 1 ? 0 : n_ % p_),
                        LCP_ + i * subarr_size, LCP_w + i * subarr_size);

            assert(is_sorted(SA_ + i * subarr_size, subarr_size + (i < p_ - 1 ? 0 : n_ % p_), LCP_ + i * subarr_size));

            if(++solved_ % 8 == 0)
                CAPS_SA_LOG(std::cerr << "\rSorted " << solved_ << " subarrays.");
        };

    parlay::parallel_for(0, p_, sort_subarr, 1);
    CAPS_SA_LOG(std::cerr << "\n");

    const auto t_e = now();
    (void) t_s; (void) t_e;
    CAPS_SA_LOG(std::cerr << "Sorted the subarrays independently. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::sort_subarrays_ext_mem()
{
    const auto t_s = now();

    std::vector<Padded<std::vector<idx_t>>> part_route_order(parlay::num_workers(), std::vector<idx_t>(p_));
    auto rng = std::default_random_engine{};
    for(auto& pd : part_route_order)
    {
        auto& v = pd.unwrap();
        std::iota(v.begin(), v.end(), 0);
        std::shuffle(v.begin(), v.end(), rng);
    }

    const auto subarr_sz = n_ / p_; // Size of each subarray to sort independently.
    const auto sort_distribute_subarr =
        [&](const idx_t p_id)
        {
            const auto w_id = parlay::worker_id();
            const auto range_beg = p_id * subarr_sz, len = subarr_sz + (p_id < p_ - 1 ? 0 : n_ % p_);
            auto& buf = w_buf[w_id].unwrap();
            buf.reserve_uninit(len);

            assert(len <= buf.SA_buf.capacity());
            auto const SA = buf.SA_buf.data(), SA_w = buf.SA_w_buf.data(), LCP = buf.LCP_buf.data(), LCP_w = buf.LCP_w_buf.data();
            for(std::size_t i = 0; i < len; ++i)
                SA[i] = SA_w[i] = range_beg + i;    // Populate a draft SA for this subproblem.

            merge_sort(SA_w, SA, len, LCP, LCP_w);
            assert(is_sorted(SA, len, LCP));

            if(++solved_ % 8 == 0)
                CAPS_SA_LOG(std::cerr << "\rSorted and partitioned " << solved_ << " subarrays.");

            // const auto pivot_off = p_id * sample_per_part_;
            // sample_pivots(SA, len, sample_per_part_, pivot_ + pivot_off);


            auto const P = buf.pivot_loc_buf.data();    // Pivot positions in this sorted subarray.

            P[0] = 0, P[p_] = len;  // Two flanking pivot indices.
            for(idx_t piv_id = 0; piv_id < p_ - 1; ++piv_id)
            {
                P[piv_id + 1] = upper_bound(SA, len, pivot_[piv_id]);
                assert(P[piv_id + 1] >= P[piv_id]);
            }

            distribute_sub_subarrays_ext_mem(part_route_order[w_id].unwrap());
        };

    solved_ = 0;
    parlay::parallel_for(0, p_, sort_distribute_subarr, 1);
    CAPS_SA_LOG(std::cerr << "\n");

    const auto t_b = now();
    parlay::parallel_for(0, p_,
    [&](const auto p_id)
    {
        subproblem_space[p_id].unwrap().close();
    }
    );

    const auto t_e = now();
    (void)t_s; (void)t_e; (void) t_b;
    CAPS_SA_LOG(std::cerr << "Sorted the subarrays independently and collated them into partitions. Time taken: " << duration(t_e - t_s) << " seconds.\n");
    CAPS_SA_LOG(std::cerr << "Closing the buckets took: " << duration(t_e - t_b) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::sample_pivots(const idx_t* const X, const idx_t n, const idx_t m, idx_t* const P)
{
    assert(m <= n);

    std::sample(X, X + n, P, m, std::mt19937(std::random_device()()));
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::select_pivots()
{
    const auto t_s = now();

    !ext_mem_ctr_ ? collect_samples() : collect_samples_ext_mem();
    select_pivots_off_samples();

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Selected the global pivots. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::collect_samples()
{
    const auto subarr_size = n_ / p_;   // Size of each sorted subarray.

    parlay::parallel_for(0, p_, [&](const idx_t p_id)
    {
        const auto samples = pivot_ + p_id * sample_per_part_;
        sample_pivots(  SA_ + p_id * subarr_size, subarr_size + (p_id < p_ - 1 ? 0 : n_ % p_),
                        sample_per_part_, samples);

        assert(is_sorted(pivot_ + p_id * sample_per_part_, sample_per_part_));
    }
    );
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::collect_samples_ext_mem()
{
    const auto subarr_sz = n_ / p_; // Size of each sorted subarray.
    std::vector<idx_t> candidates_off(subarr_sz + n_ % p_);
    std::iota(candidates_off.begin(), candidates_off.end(), idx_t(0));

    parlay::parallel_for(0, p_, [&](const idx_t p_id)
    {
        const auto samples = pivot_ + p_id * sample_per_part_;
        sample_pivots(  candidates_off.data(), subarr_sz + (p_id < p_ - 1? 0 : n_ % p_),
                        sample_per_part_, samples);

        const auto samples_off = p_id * subarr_sz;
        std::for_each(pivot_ + p_id * sample_per_part_, pivot_ + (p_id + 1) * sample_per_part_, [&](auto& s){ s += samples_off; });
    }
    );
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::select_pivots_off_samples()
{
    const auto sample_count = p_ * sample_per_part_;    // Total number of samples to select pivots from.
    idx_t* const pivot_w = allocate<idx_t>(sample_count);   // Working space to sample pivots.
    auto const temp_1 = allocate<idx_t>(sample_count), temp_2 = allocate<idx_t>(sample_count);

    std::memcpy(pivot_w, pivot_, sample_count * sizeof(idx_t));
    merge_sort(pivot_, pivot_w, sample_count, temp_1, temp_2);
    assert(is_sorted(pivot_w, sample_count, temp_1));

    const auto gap = sample_count / (p_);   // Distance-gap between pivots.
    const auto quant_err = sample_count % p_;   // Numerator of the quantization error per partition.
    idx_t unaccounted = 0;  // Running total number of samples not properly accounted for due to quantization of partition sizes.
    std::size_t idx = 0;    // Index of the next pivot to choose from the samples.
    for(idx_t i = 0; i < p_ - 1; ++i)
    {
        idx += gap;
        unaccounted += quant_err;
        if(unaccounted >= p_)   // Fix the cumulative error.
            idx++,
            unaccounted %= p_;

        pivot_[i] = pivot_w[idx - 1];
    }

    deallocate(pivot_w), deallocate(temp_1), deallocate(temp_2);
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::locate_pivots(idx_t* const P) const
{
    const auto t_s = now();

    const auto subarr_size = n_ / p_;   // Size of each independent sorted subarray.

    const auto locate =
        [&](const idx_t i)
        {
            const auto X_i = SA_ + i * subarr_size; // The i'th subarray.
            const auto P_i = P + i * (p_ + 1);  // Pivot locations in `X_i` are to be placed in `P_i`.

            const auto tot_subarr_size = subarr_size + (i < p_ - 1 ? 0 : n_ % p_);
            P_i[0] = 0, P_i[p_] = tot_subarr_size; // The two flanking pivot indices.

            for(idx_t j = 0; j < p_ - 1; ++j) // TODO: try parallelizing this loop too; observe performance diff.
                P_i[j + 1] = upper_bound(X_i, tot_subarr_size, pivot_[j]);  // TODO: can shrink the search-range for successive searches.
        };

    parlay::parallel_for(0, p_, locate, 1);

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Located the pivots in each sorted subarray. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
T_idx_ Suffix_Array<T_seq_, T_idx_>::upper_bound(const idx_t* const X, const idx_t n, const idx_t p) const
{
    // Invariant: SA[l] < P < SA[r].

    const auto P_len = n_ - p;
    int64_t l = -1, r = n;  // (Exclusive-) Range of the iterations in the binary search.
    idx_t c;    // Midpoint in each iteration.
    idx_t soln = n; // Solution of the search.
    idx_t lcp_l = 0, lcp_r = 0; // LCP(P, SA[l]) and LCP(P, SA[r]).
	constexpr idx_t cutoff = 65536; // TODO: better tune and document.

    while(r - l > 1)    // Candidate matches exist.
    {
        c = (l + r) / 2;
        auto const suf = X[c];  // The suffix at the middle.
        const auto suf_len = n_ - suf;  // Length of the suffix.

        idx_t lcp_c = std::min(lcp_l, lcp_r);   // LCP(X[c], P).
        lcp_c = std::min(lcp_c, cutoff);
        auto max_lcp = std::min(std::min(suf_len, P_len), max_context); // Maximum possible LCP, i.e. length of the shorter string.
        max_lcp = std::min(max_lcp, cutoff);
        lcp_c += LCP(suf + lcp_c, p + lcp_c, max_lcp - lcp_c);  // Skip an informed number of character comparisons.

        if(lcp_c == max_lcp)    // One is a prefix of the other, or they align at least up-to the context- or the cutoff-length.
        {
            if(lcp_c == P_len)  // P is a prefix of the suffix.
            {
                if(P_len == suf_len)    // The query is the suffix itself, i.e. P = X[c]
                    return c + 1;
                else    // P < X[c]
                    r = c, lcp_r = lcp_c, soln = c;
            }
            else    // The suffix is a prefix of the query, so X[c] < P; technically impossible if the text terminates with $.
                    // Or, their relevant prefixes align; moving to the right, as searching for the upper-bound.
                l = c, lcp_l = lcp_c;
        }
        else    // They mismatch within their relevant prefixes.
            if(is_lesser(suf + lcp_c, p + lcp_c))   // X[c] < P
                l = c, lcp_l = lcp_c;
            else    // P < X[c]
                r = c, lcp_r = lcp_c, soln = c;
    }


    return soln;
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::partition_sub_subarrays(const idx_t* const P)
{
    const auto t_s = now();

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
    prefix_sum(part_size_scan_, p_);
    assert(part_size_scan_[p_] == n_);


    // Collate the sorted sub-subarrays to appropriate partitions.
    const idx_t subarr_size = n_ / p_;
    const auto collate =    // Collates the `j`'th sub-subarray from each sorted subarray to partition `j`.
        [&](const idx_t j)
        {
            auto const Y_j = SA_w + part_size_scan_[j]; // Memory-base for partition `j`.
            auto const LCP_Y_j = LCP_w + part_size_scan_[j];    // Memory-base for LCPs of partition `j`.
            auto const sub_subarr_off = part_ruler_ + j * (p_ + 1); // Offset of the sorted sub-subarrays in `Y_j`.
            idx_t curr_idx = 0; // Current index into `Y_j`.

            for(idx_t i = 0; i < p_; ++i)   // Subarray `i`.
            {
                const auto X_i = SA_ + i * subarr_size; // `i`'th sorted subarray.
                const auto LCP_X_i = LCP_ + i * subarr_size;    // LCP array of `X_i`.
                const auto P_i = P + i * (p_ + 1);  // Pivot collection of subarray `i`.

                sub_subarr_off[i] = curr_idx;
                const auto sub_subarr_size = P_i[j + 1] - P_i[j];   // Size of the `j`'th sub-subarray of subarray `i`.
                if(sub_subarr_size == 0)
                    continue;

                std::memcpy(Y_j + sub_subarr_off[i], X_i + P_i[j], sub_subarr_size * sizeof(idx_t));
                std::memcpy(LCP_Y_j + sub_subarr_off[i], LCP_X_i + P_i[j], sub_subarr_size * sizeof(idx_t));
                LCP_Y_j[sub_subarr_off[i]] = 0;
                curr_idx += sub_subarr_size;
            }

            sub_subarr_off[p_] = curr_idx;
            assert(curr_idx == part_size_scan_[j + 1] - part_size_scan_[j]);
        };

    parlay::parallel_for(0, p_, collate, 1);

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Collated the sorted sub-subarrays into partitions. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::distribute_sub_subarrays_ext_mem(const std::vector<idx_t>& route_order)
{
    const auto w_id = parlay::worker_id();
    auto& buf = w_buf[w_id].unwrap();
    auto const SA = buf.SA_buf.data(), LCP = buf.LCP_buf.data(), P = buf.pivot_loc_buf.data();


    // Different sequences of parts-copying (dispersion) by different workers to minimize lock-contention.
    for(const auto part_id : route_order)   // TODO: why not in parallel?
    {
        const auto sub_subarr_sz = P[part_id + 1] - P[part_id];
        if(sub_subarr_sz > 0)
            LCP[P[part_id]] = 0;    // An independent sorted segment's first LCP-value is 0.

        assert(is_sorted(SA + P[part_id], sub_subarr_sz, LCP + P[part_id]));

        lock[part_id].unwrap().lock();

        auto& b = subproblem_space[part_id].unwrap();
        b.add(SA + P[part_id], LCP + P[part_id], sub_subarr_sz);

        lock[part_id].unwrap().unlock();
    }
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::merge_sub_subarrays()
{
    const auto t_s = now();

    const auto dup =    // Duplicates the `j`'th partition.
        [&](const idx_t j)
        {
            const auto part_size = part_size_scan_[j + 1] - part_size_scan_[j];
            std::memcpy(SA_ + part_size_scan_[j], SA_w + part_size_scan_[j], part_size * sizeof(idx_t));
            std::memcpy(LCP_ + part_size_scan_[j], LCP_w + part_size_scan_[j], part_size * sizeof(idx_t));
        };

    parlay::parallel_for(0, p_, dup, 1);    // Fulfill `sort_partition`'s precondition.


    std::atomic_uint64_t solved_ = 0;   // Progress tracker—number of subproblems solved in some step.
    const auto sort_part =
        [&](const idx_t j)
        {
            const auto part_off = part_size_scan_[j];   // Offset of the partition in the partitions' flat collection.
            auto const X_j = SA_w + part_off;   // Memory-base for partition `j`.
            auto const Y_j = SA_ + part_off;    // Location to sort partition `j`.
            auto const LCP_X_j = LCP_w + part_off;  // Memory-base for the LCP-arrays of partition `j`.
            auto const LCP_Y_j = LCP_ + part_off;   // LCP array of `Y_j`.
            auto const sub_subarr_off = part_ruler_ + j * (p_ + 1); // Indices of the sorted subarrays in `X_i`.

            sort_partition(X_j, Y_j, p_, sub_subarr_off, LCP_X_j, LCP_Y_j);

            if(++solved_ % 8 == 0)
               CAPS_SA_LOG(std::cerr << "\rMerged " << solved_ << " partitions.");
        };

    parlay::parallel_for(0, p_, sort_part, 1);  // Merge the sorted subarrays in each partitions.
    CAPS_SA_LOG(std::cerr << "\n");

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Merged the sorted subarrays in each partition. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::merge_sub_subarrays_ext_mem()
{
    const auto t_s = now();

    std::vector<Padded<idx_t*>> sub_subarr_idx_buf(parlay::num_workers());  // Buffer for the sorted sub-subarray sizes in each partition for each worker.
    std::for_each(sub_subarr_idx_buf.begin(), sub_subarr_idx_buf.end(), [this](auto& buf){ buf.unwrap() = allocate<idx_t>(p_ + 1); });

    const auto merge_sub_subarray =
        [&](const idx_t p_id)
        {
            const auto w_id = parlay::worker_id();
            auto& buf = w_buf[w_id].unwrap();
            auto sub_subarr_idx = sub_subarr_idx_buf[w_id].unwrap();

            auto& b = subproblem_space[p_id].unwrap();
            auto &SA_b = b.SA_bucket, &LCP_b = b.LCP_bucket;

            const idx_t part_sz = SA_b.size();
            buf.reserve_uninit(part_sz);

            auto &SA = buf.SA_buf, &LCP = buf.LCP_buf, &SA_w = buf.SA_w_buf, &LCP_w = buf.LCP_w_buf;

            // Load buckets and fulfill `sort_partition`'s precondition.
            b.load(SA.data(), LCP.data(), sub_subarr_idx + 1);
            std::memcpy(SA_w.data(), SA.data(), part_sz * sizeof(idx_t));
            std::memcpy(LCP_w.data(), LCP.data(), part_sz * sizeof(idx_t));

            sub_subarr_idx[0] = 0;
            for(idx_t i = 1; i <= p_; ++i)  // Convert sub-subarray sizes to sub-subarray indices.
                sub_subarr_idx[i] += sub_subarr_idx[i - 1];

            assert(sub_subarr_idx[p_] == part_sz);

            for(idx_t i = 0; i < p_; ++i)
                assert(is_sorted(SA.data() + sub_subarr_idx[i], sub_subarr_idx[i + 1] - sub_subarr_idx[i], LCP.data() + sub_subarr_idx[i]));

            sort_partition(SA_w.data(), SA.data(), p_, sub_subarr_idx, LCP_w.data(), LCP.data());
            assert(is_sorted(SA.data(), part_sz, LCP.data()));

            SA_b.rewrite(SA.data(), part_sz);
            if(op_lcp)  // TODO: output a json file and note this.
                LCP_b.rewrite(LCP.data(), part_sz); // TODO: note that `LCP[0] = 0`, which needs to be updated when concatenating the partitions afterwards.

            if(++solved_ % 8 == 0)
                CAPS_SA_LOG(std::cerr << "\rMerged " << solved_ << " partitions.");
        };

    solved_ = 0;
    parlay::parallel_for(0, p_, merge_sub_subarray, 1);
    CAPS_SA_LOG(std::cerr << "\n");


    std::for_each(sub_subarr_idx_buf.cbegin(), sub_subarr_idx_buf.cend(), [](auto& buf){ deallocate(buf.unwrap()); });

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Merged the sorted subarrays in each partition. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::sort_partition(idx_t* const X, idx_t* const Y, const idx_t n, const idx_t* const S, idx_t* const LCP_x, idx_t* const LCP_y)
{
    if(n == 1)
        return;

    const auto m = n / 2;
    const auto flat_count_l = S[m] - S[0];
    const auto flat_count_r = S[n] - S[m];

    const auto f = [&](){ sort_partition(Y, X, m, S, LCP_y, LCP_x); };
    const auto g = [&](){ sort_partition(Y + flat_count_l, X + flat_count_l, n - m, S + m, LCP_y + flat_count_l, LCP_x + flat_count_l); };

    (flat_count_l < nested_par_grain_size || flat_count_r < nested_par_grain_size) ?
        (f(), g()) : parlay::par_do(f, g, ext_mem_ctr_);
    merge(X, flat_count_l, X + flat_count_l, flat_count_r, LCP_x, LCP_x + flat_count_l, Y, LCP_y);
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::compute_partition_boundary_lcp()
{
    const auto t_s = now();

    const auto compute_boundary_lcp =
        [&](const idx_t j)
        {
            const auto part_idx = part_size_scan_[j];
            LCP_[part_idx] = LCP(SA_[part_idx - 1], SA_[part_idx], n_ - std::max(SA_[part_idx - 1], SA_[part_idx]));
        };

    parlay::parallel_for(1, p_, compute_boundary_lcp, 1);

    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Computed the LCPs at the partition boundaries. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}

template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::remove_extmem_partitions()
{
    if(!ext_mem_ctr_) {
        CAPS_SA_LOG(std::cerr << "The remove_extmem_partitions() member has no effect when not using external memory construction. Please ensure you intended to call this member.\n");
        return;
    } 
    
    if(!constructed) {
        CAPS_SA_LOG(std::cerr << "The remove_extmem_partitions() member cannot be called on a Suffix_Array object that has not yet been constructed. Please ensure you intended to call this member here.\n");
        return;
    }
    if(!removed_extmem_partitions) {
        const auto remove_bucket = [this](const idx_t p_id){ 
            subproblem_space[p_id].unwrap().remove_SA(); 
            if(!removed_lcp_partitions) {
                subproblem_space[p_id].unwrap().remove_LCP(); 
            }
        };
        parlay::parallel_for(0, p_, remove_bucket, 1);
        removed_extmem_partitions = true;
    } else {
        CAPS_SA_LOG(std::cerr << "External memory partitions have already been removed. Please check why you may be calling this member more than once.\n");
    }
}

template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::clean_up()
{
    const auto t_s = now();

    if(!ext_mem_ctr_)
    {
        deallocate(SA_w), deallocate(LCP_w);

        deallocate(pivot_);
        deallocate(part_size_scan_);
        deallocate(part_ruler_);
    }
    else
    {
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            w_buf[w_id].unwrap().free();

        deallocate(pivot_);
        const auto clean_bucket = [this](const idx_t p_id){ 
            if(!op_lcp) {
                subproblem_space[p_id].unwrap().remove_LCP();
            }
            subproblem_space[p_id].unwrap().clear(); 
        };
        parlay::parallel_for(0, p_, clean_bucket, 1);
        removed_lcp_partitions = !op_lcp;
    }
    
    constructed = true;
    const auto t_e = now();
    (void)t_s; (void)t_e;
    CAPS_SA_LOG(std::cerr << "Released the temporary data structures. Time taken: " << duration(t_e - t_s) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::print_stats() const
{
    std::vector<uint64_t> p_sz(p_); // Partition sizes.
    for(std::size_t p_id = 0; p_id < p_; ++p_id)
        p_sz[p_id] = !ext_mem_ctr_ ?
                        part_size_scan_[p_id + 1] - part_size_scan_[p_id] :
                        subproblem_space[p_id].unwrap().SA_bucket.size();

    CAPS_SA_LOG(std::cerr << "Bucket stats: " << "\n");

    const auto sum =  std::accumulate(p_sz.cbegin(), p_sz.cend(), uint64_t(0));
    const auto mean = static_cast<double>(sum) / p_;
    double var = 0;
    std::for_each(p_sz.cbegin(), p_sz.cend(), [&](const auto sz){ var += (sz - mean) * (sz - mean); });
    var /= p_;
    const auto sd = std::sqrt(var);

    (void) sd;
    CAPS_SA_LOG(std::cerr << "\t Sum size:  " << sum << "\n");
    CAPS_SA_LOG(std::cerr << "\t Max size:  " << *std::max_element(p_sz.cbegin(), p_sz.cend()) << "\n");
    CAPS_SA_LOG(std::cerr << "\t Min size:  " << *std::min_element(p_sz.cbegin(), p_sz.cend()) << "\n");
    CAPS_SA_LOG(std::cerr << "\t Mean size: " << mean << "\n");
    CAPS_SA_LOG(std::cerr << "\t SD(size):  " << sd << "\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::construct()
{
    const auto t_start = now();

    initialize();

    permute();

    // merge_sort(SA_w, SA_, n_, LCP_, LCP_w);  // Monolithic construction.

    sort_subarrays();

    select_pivots();

    idx_t* const P = allocate<idx_t>(p_ * (p_ + 1));  // Collection of pivot locations in the subarrays.
    locate_pivots(P);
    partition_sub_subarrays(P);
    deallocate(P);

    print_stats();

    merge_sub_subarrays();

    compute_partition_boundary_lcp();

    clean_up();

    constructed = true;
    const auto t_end = now();
    (void)t_start; (void)t_end;
    CAPS_SA_LOG(std::cerr << "Constructed the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::construct_ext_mem()
{
    const auto t_start = now();

    initialize();

    select_pivots();

    sort_subarrays_ext_mem();
    // select_pivots_off_samples();

    print_stats();

    merge_sub_subarrays_ext_mem();

    clean_up();

    constructed = true;
    const auto t_end = now();
    (void)t_start; (void)t_end;
    CAPS_SA_LOG(std::cerr << "Constructed the suffix array and the LCP array. Time taken: " << duration(t_end - t_start) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
const T_idx_* Suffix_Array<T_seq_, T_idx_>::SA()
{
    if(SA_ == nullptr)
    {
        assert(ext_mem_ctr_);

        idx_t* const pref_sum = allocate<idx_t>(p_ + 1);
        parlay::parallel_for(0, p_,
            [&](const std::size_t idx){ pref_sum[idx] = subproblem_space[idx].unwrap().SA_bucket.size(); });

        prefix_sum(pref_sum, p_);


        SA_ = allocate<idx_t>(n_);
        parlay::parallel_for(0, p_,
            [&](const std::size_t idx)
            {
                const auto& SA_bucket = subproblem_space[idx].unwrap().SA_bucket;
                const auto read_elems = SA_bucket.load(SA_ + pref_sum[idx]);
                assert(read_elems == SA_bucket.size());
                (void)read_elems;
            });

        deallocate(pref_sum);
    }

    return SA_;
}


template <typename T_seq_, typename T_idx_>
void Suffix_Array<T_seq_, T_idx_>::dump(std::ofstream& output) const 
{
    const auto t_start = now();

    const std::size_t n = n_;
    output.write(reinterpret_cast<const char*>(&n), sizeof(std::size_t));

    if(!ext_mem_ctr_)
    {
        output.write(reinterpret_cast<const char*>(SA_), n_ * sizeof(idx_t));
        if(op_lcp)
            output.write(reinterpret_cast<const char*>(LCP_), n_ * sizeof(idx_t));
    }
    else
    {
        for(idx_t p_id = 0; p_id < p_; ++p_id)
        {
            std::ifstream input(SA_bucket_file_path(p_id));
            assert(input.peek() != EOF);
            output << input.rdbuf();
            input.close();
        }

        if(op_lcp)
            for(idx_t p_id = 0; p_id < p_; ++p_id)
            {
                std::ifstream input(LCP_bucket_file_path(p_id));
                assert(input.peek() != EOF);
                output << input.rdbuf();
                input.close();
            }

        output.close();
    }

    const auto t_end = now();
    (void)t_start; (void)t_end;
    CAPS_SA_LOG(std::cerr << "Dumped the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n");
}


template <typename T_seq_, typename T_idx_>
template <typename T_arr_>
void Suffix_Array<T_seq_, T_idx_>::prefix_sum(T_arr_* A, const T_idx_ n)
{
    T_arr_ curr_sum = 0;
    for(idx_t i = 0; i < n; ++i)
    {
        const auto curr_val = A[i];

        A[i] = curr_sum;
        curr_sum += curr_val;
    }

    A[n] = curr_sum;
}


template <typename T_seq_, typename T_idx_>
bool Suffix_Array<T_seq_, T_idx_>::is_smaller(const idx_t x, const idx_t y, idx_t& lcp) const
{
    assert(x < n_ && y < n_);

    const auto l = std::min(n_ - x, n_ - y);    // TODO: add context-bound.

    // Using byte-by-byte variant for correctness-check.
    if constexpr(std::is_same_v<T_seq_, Genomic_Text>)
        lcp = T_->lcp_unvectorized(x, y, l);
    else
        lcp = lcp_unvectorized(reinterpret_cast<const char*>(T_ + x), reinterpret_cast<const char*>(T_ + y), l * sizeof(T_seq_)) / sizeof(T_seq_);  // Using byte-by-byte variant for correctness-check.

    if(lcp == l)    // Shorter suffix is a prefix of the longer one;
                    // possible in strings w/o designated delimiters, and in context-bounded SAs.
    {
        return x > y;   // Shorter suffixes are considered lexicographically smaller in such cases.
    }

    return is_lesser(x + lcp, y + lcp);
}


template <typename T_seq_, typename T_idx_>
bool Suffix_Array<T_seq_, T_idx_>::is_sorted(const idx_t* const X, const idx_t n, const idx_t* const LCP_X) const
{
    if(LCP_X && n > 0 && LCP_X[0] != 0)
        return false;

    idx_t lcp;
    for(idx_t i = 1; i < n; ++i)
        if(!is_smaller(X[i - 1], X[i], lcp))
            return false;
        else if(LCP_X && lcp != LCP_X[i])
            return false;

    return true;
}


template <typename T_seq_, typename T_idx_>
bool Suffix_Array<T_seq_, T_idx_>::is_correct()
{
    std::vector<Padded<bool>> result(p_, true);    // Correctness result from each worker.

    if(SA_ == nullptr)  // External-memory scenario.
    {
        assert(ext_mem_ctr_);

        std::vector<Padded<idx_t*>> SA_buf_w(parlay::num_workers());    // Buffer space for each worker, to load SA parts.
        std::vector<Padded<idx_t*>> LCP_buf_w(parlay::num_workers());   // Buffer space for each worker, to load LCP-array parts.

        std::vector<Padded<idx_t>> s_first(p_); // First suffix in each part.
        std::vector<Padded<idx_t>> s_last(p_);  // Last suffix in each part.

        std::size_t max_b_size = 0; // Maximum bucket-size.
        std::for_each(subproblem_space.cbegin(), subproblem_space.cend(), [&](const auto& b){ max_b_size = std::max(max_b_size, b.unwrap().SA_bucket.size()); });

        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            SA_buf_w[w_id] = allocate<idx_t>(max_b_size),
            LCP_buf_w[w_id] = allocate<idx_t>(max_b_size);


        // Method to check if a part is sorted.
        const auto check_sorted = [&](const std::size_t p_id)
        {
            const auto& b = subproblem_space[p_id].unwrap();
            const auto& SA_b = b.SA_bucket;
            const auto& LCP_b = b.LCP_bucket;
            const auto sz = SA_b.size();
            assert(sz == LCP_b.size());
            auto SA_buf = SA_buf_w[parlay::worker_id()].unwrap(), LCP_buf = LCP_buf_w[parlay::worker_id()].unwrap();

            SA_b.load(SA_buf);
            LCP_b.load(LCP_buf);

            if(!is_sorted(SA_buf, sz))//, LCP_buf)) // TODO: fix this—LCP buckets may not have updated data.
                result[p_id].unwrap() = false;

            s_first[p_id].unwrap() = SA_buf[0], s_last[p_id].unwrap() = SA_buf[sz - 1];
        };

        parlay::parallel_for(0, p_, check_sorted, 1);


        // Method to check if the parts are in sorted order.
        const auto check_part_order = [&](const std::size_t p_id)
        {
            if(p_id == 0)
                return;

            idx_t lcp;
            if(!is_smaller(s_last[p_id - 1].unwrap(), s_first[p_id].unwrap(), lcp))
                result[p_id].unwrap() = false;

            // Note: we don't have the inter-part (i.e. boundary) LCP-information yet—so no use of `lcp`.
        };

        parlay::parallel_for(0, p_, check_part_order, 1);


        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            deallocate(SA_buf_w[w_id].unwrap()), deallocate(LCP_buf_w[w_id].unwrap());
    }
    else
    {
        if(LCP_ && LCP_[0] != 0)
            return false;

        const auto range_sz_per_p = n_ / p_;
        const auto is_sorted = [&](const std::size_t p_id)
        {
            const auto beg = (p_id > 0 ? range_sz_per_p * p_id : 1);
            const auto end = (p_id < p_ - 1 ? (p_id + 1) * range_sz_per_p : n_);

            idx_t lcp;
            for(idx_t i = beg; i < end; ++i)
                if(!is_smaller(SA_[i - 1], SA_[i], lcp))
                    result[p_id].unwrap() = false;
                else if(LCP_ && lcp != LCP_[i])
                    result[p_id].unwrap() = false;
        };

        parlay::parallel_for(0, p_, is_sorted, 1);

        // The parts' order is implicitly checked by comparing `i - 1` and `i`.
        // No additional inter-part comparisons required.
    }


    bool correct = true;
    std::for_each(result.cbegin(), result.cend(), [&correct](const auto& r){ correct = correct && r.unwrap(); });
    return correct;
}

}



// Template instantiations for the required instances.
template class CaPS_SA::Suffix_Array<char, uint32_t>;
template class CaPS_SA::Suffix_Array<char, uint64_t>;
template class CaPS_SA::Suffix_Array<uint32_t, uint32_t>;
template class CaPS_SA::Suffix_Array<uint32_t, uint64_t>;
template class CaPS_SA::Suffix_Array<uint64_t, uint32_t>;
template class CaPS_SA::Suffix_Array<uint64_t, uint64_t>;
template class CaPS_SA::Suffix_Array<CaPS_SA::Genomic_Text, uint32_t>;
template class CaPS_SA::Suffix_Array<CaPS_SA::Genomic_Text, uint64_t>;
