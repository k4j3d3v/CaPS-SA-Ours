
#ifndef THEMIS_SUFFIX_ARRAY_HPP
#define THEMIS_SUFFIX_ARRAY_HPP



#include "Ext_Mem_Bucket.hpp"
#include "Spin_Lock.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <immintrin.h>

// =============================================================================

namespace CaPS_SA
{

// The Suffix Array (SA) and the Longest Common Prefix (LCP) array constructor
// class for some given sequence.
template <typename T_idx_>
class Suffix_Array
{
private:

    typedef T_idx_ idx_t;   // Integer-type for indexing the input text.
    typedef Padded<idx_t*> w_local_buf_t;   // Type of worker-local buffers.

    const char* const T_;   // The input text.
    const idx_t n_; // Length of the input text.

    const idx_t p_; // Count of subproblems used in construction.

    const idx_t per_worker_in_mem_elem; // Maximum number of elements for each worker to keep in memory in external-memory setting.

    idx_t* SA_; // The SA in internal memory.
    idx_t* LCP_;    // The LCP-array in internal memory.

    idx_t* SA_w;    // Working space for the SA construction in internal memory.
    idx_t* LCP_w;   // Working space for the LCP-array construction in internal memory.

    // In-memory space for a worker to operate on a subproblem, in external-memory setting.
    struct Worker_Mem;

    // External-memory space corresponding to a subproblem.
    struct Subproblem_Ext_Mem;

    std::vector<Padded<Worker_Mem>> w_buf;  // In-memory buffer for each worker.
    std::vector<Padded<Subproblem_Ext_Mem>> subproblem_space;   // External-memory space for each subproblem.

    std::vector<Padded<Spin_Lock>> lock;    // Lock for each external-memory bucket. TODO: use the concurrent ext-mem bucket impl, circumventing a whole-bucket lock.

    const bool ext_mem_ctr_;    // Whether to construct using external-memory or not.

    const std::string ext_mem_path; // Path-prefix to external-memory buckets.

    const idx_t max_context;    // Maximum prefix-context length for comparing suffixes.

    const idx_t sample_per_part_; // Number of pivots to sample per subarray.
    idx_t* pivot_;  // Pivots for the global SA.
    idx_t* part_size_scan_; // Inclusive scan (prefix sum) of the sizes of the pivoted final partitions containing appropriate sorted sub-subarrays.
    idx_t* part_ruler_; // "Ruler" for the partitions—contains the indices of each sub-subarray in each partition.

    std::atomic_uint64_t solved_;   // Progress tracker—number of subproblems solved.

    const bool op_lcp;  // Whether to output the LCP-array.

    static constexpr idx_t default_subproblem_count = 8192; // Default subproblem-count to use in construction.
    static constexpr idx_t nested_par_grain_size = (1lu << 13); // Granularity for nested parallelism to kick in.


    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`.
    static idx_t lcp(const char* x, const char* y, idx_t min_len);

    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`. `N x 32` bytes of prefix comparisons are
    // loop-unrolled.
    template <std::size_t N = 8>
    static idx_t LCP(const char* x, const char* y, idx_t min_len);

    // Returns the LCP length of the `32 x N`-bytes prefix of `x` and `y`.
    template <std::size_t N> static idx_t LCP_unrolled(const char* x, const char* y);

    // Merges the sorted collections of suffixes, `X` and `Y`, with lengths
    // `len_x` and `len_y` and LCP-arrays `LCP_x` and `LCP_y` respectively, into
    // `Z`. Also constructs `Z`'s LCP-array in `LCP_z`.
    void merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const;

    // Merge-sorts the suffix collection `X` of length `n` into `Y`. Also
    // constructs the LCP-array of `X` in `LCP`, using `W` as working space.
    // A necessary precondition is that `Y` must be equal to `X`.  `X` may
    // not remain the same after the sort.
    void merge_sort(idx_t* X, idx_t* Y, idx_t n, idx_t* LCP, idx_t* W) const;

    // Initializes internal data structures for the construction algorithm.
    void initialize();

    // Populates a draft SA with some permutation of `[0, len)`.
    void permute();

    // Sorts uniform-sized subarrays independently.
    void sort_subarrays();

    // Sorts uniform-sized subarrays independently in external-memory setting.
    void sort_subarrays_ext_mem();

    // Samples `m` pivots from the sorted suffix collection `X` of size `n`
    // into `P`.
    static void sample_pivots(const idx_t* X, idx_t n, idx_t m, idx_t* P);

    // Samples `m` pivots from the suffix collection range `[r_beg, r_beg + n)`
    // into `P`.
    static void sample_pivots(const idx_t r_beg, idx_t n, idx_t m, idx_t* P);

    // Selects pivots for parallel merging of the sorted subarrays.
    void select_pivots();

    // Collects samples from each sorted subarray.
    void collect_samples();

    // Collects samples from each subarray, for the external-memory setting of
    // the algorithm.
    void collect_samples_ext_mem();

    // Selects the final pivot set from the sampled suffixes.
    void select_pivots_off_samples();

    // Locates the positions (upper-bounds) of the selected pivots in the sorted
    // subarrays and flattens them in `P`. Besides these pivots, two flanking
    // pivots, `0` and `|X_i|`, for each subarray `X_i` are also placed.
    void locate_pivots(idx_t* P) const;

    // Returns the first index `idx` into the sorted suffix collection `X` of
    // length `n` such that `X[idx]` is strictly greater than the query pattern
    // `P` of length `P_len`.
    idx_t upper_bound(const idx_t* X, idx_t n, const char* P, idx_t P_len) const;

    // Collates the sub-subarrays delineated by the pivot locations in each
    // sorted subarray, present in `P`, into appropriate partitions.
    void partition_sub_subarrays(const idx_t* P);

    // Distributes the sub-subarrays delineated by the pivot locations in the
    // current sorted subarray being processed, into appropriate partitions (in
    // external-memory). The sub-subarrays are distributed as per the order in
    // `route_order`.
    void distribute_sub_subarrays_ext_mem(const std::vector<idx_t>& route_order);

    // Merges the sorted sub-subarrays laid flat together in each partition.
    void merge_sub_subarrays();

    // Merges the sorted sub-subarrays in each external-memory partition.
    void merge_sub_subarrays_ext_mem();

    // Computes the LCPs at the partition boundaries, specifically at the
    // starting index of each partition in their flat collection.
    void compute_partition_boundary_lcp();

    // Merge-sorts the collection `X` that contains `n` sorted arrays of
    // suffixes laid flat together, into `Y`. `S` contains the delimiter indices
    // of the `n` arrays in `X`. The LCP-array of sorted `X` is constructed in
    // `LCP_y`; `LCP_x` contains the LCP-arrays of the `n` arrays of `X`.
    // A necessary precondition is that `Y` must be equal to `X`, and `LCP_y` to
    // `LCP_x`. `X` and `LCP_x` may not remain the same after the sort.
    void sort_partition(idx_t* X, idx_t* Y, idx_t n, const idx_t* S, idx_t* LCP_x, idx_t* LCP_y);

    // Cleans up after the construction algorithm.
    void clean_up();

    // Cleans up after the external-memory construction algorithm.
    void clean_up_ext_mem();

    const std::string SA_bucket_file_path(const idx_t p_id) const { return ext_mem_path + "_SA_" + std::to_string(p_id); }
    const std::string LCP_bucket_file_path(const idx_t p_id) const { return ext_mem_path + "_LCP_" + std::to_string(p_id); }
    const std::string sz_bucket_file_path(const idx_t p_id) const { return ext_mem_path + "_sz_" + std::to_string(p_id); }

    // Computes in-place prefix-sum in the array `A` for the first `n` elements.
    // `A` must have size at least `n + 1`—one extra entry is required to hold
    // the total sum.
    template <typename T_>
    static void prefix_sum(T_* A, idx_t n);

    // Returns `true` iff the suffix `x` is lexicographically smaller than the
    // suffix `y`. Their LCP is computed in `lcp`.
    bool is_smaller(idx_t x, idx_t y, idx_t& lcp) const;

    // Returns true iff `X` is a valid (partial) SA with size `n`. An optional
    // LCP-array `LCP_X` can be provided to check its correctness too.
    bool is_sorted(const idx_t* X, idx_t n, const idx_t* LCP_X = nullptr) const;

public:

    // Constructs a suffix array object for the input text `T` of size
    // `n`. External-memory is used for construction if `ext_mem` is specified.
    // For the external-memory setting, the path-prefix `ext_mem_path` is used.
    // Optionally, the number of subproblems to decompose the original
    // construction problem into can be provided with `subproblem_count`, and
    // the maximum prefix-context length for the suffixes can be bounded by
    // `max_context`.
    Suffix_Array(const char* T, idx_t n, bool ext_mem = false, const std::string& ext_mem_path = ".", idx_t subproblem_count = 0, idx_t max_context = 0);

    Suffix_Array(const Suffix_Array&) = delete;
    Suffix_Array& operator=(const Suffix_Array&) = delete;
    Suffix_Array(Suffix_Array&&) = delete;
    Suffix_Array& operator=(Suffix_Array&&) = delete;

    ~Suffix_Array();

    // Returns the text.
    const char* T() const { return T_; }

    // Returns the length of the text.
    idx_t n() const { return n_; }

    // Returns the SA.
    const idx_t* SA();

    // Returns the LCP-array.
    const idx_t* LCP() const { return LCP_; }

    // Constructs the SA and the LCP-array in internal memory.
    void construct();

    // Constructs the SA and the LCP-array using external memory.
    void construct_ext_mem();

    // Dumps the SA and the LCP-array into the stream `output`.
    void dump(std::ofstream& output) const;

    // Returns `true` iff the SA (and the LCP-array, if retained) are
    // correct.
    bool is_correct();
};


template <typename T_idx_>
struct Suffix_Array<T_idx_>::Worker_Mem
{
    // TODO: use custom `Buffer` impl.

    idx_t* SA_buf;  // Memory buffer for SA elements.
    idx_t* LCP_buf; // Memory buffer for LCP-array elements

    idx_t* SA_w_buf;    // Working space for the SA construction.
    idx_t* LCP_w_buf;   // Working space for the LCP-array construction.

    idx_t* pivot_loc_buf;   // Buffers to store pivot locations in the sorted subarray.


    Worker_Mem(const Suffix_Array<T_idx_>& SA):
          SA_buf(allocate<idx_t>(SA.per_worker_in_mem_elem))
        , LCP_buf(allocate<idx_t>(SA.per_worker_in_mem_elem))
        , SA_w_buf(allocate<idx_t>(SA.per_worker_in_mem_elem))
        , LCP_w_buf(allocate<idx_t>(SA.per_worker_in_mem_elem))
        , pivot_loc_buf(allocate<idx_t>(SA.p_ + 2))
    {}

    Worker_Mem(Worker_Mem&& rhs):
          SA_buf(rhs.SA_buf)
        , LCP_buf(rhs.LCP_buf)
        , SA_w_buf(rhs.SA_w_buf)
        , LCP_w_buf(rhs.LCP_w_buf)
        , pivot_loc_buf(rhs.pivot_loc_buf)
    {
        rhs.SA_buf = rhs.LCP_buf = rhs.SA_w_buf = rhs.LCP_w_buf = rhs.pivot_loc_buf = nullptr;
    }

    Worker_Mem(const Worker_Mem&) = delete;
    Worker_Mem& operator=(const Worker_Mem&) = delete;
    Worker_Mem&& operator=(Worker_Mem&&) = delete;

    void free()
    {
        deallocate(SA_buf), deallocate(LCP_buf), deallocate(SA_w_buf), deallocate(LCP_w_buf), deallocate(pivot_loc_buf);
        SA_buf = LCP_buf = SA_w_buf = LCP_w_buf = pivot_loc_buf = nullptr;
    }

    ~Worker_Mem()
    {
        deallocate(SA_buf), deallocate(LCP_buf), deallocate(SA_w_buf), deallocate(LCP_w_buf), deallocate(pivot_loc_buf);
    }
};


template <typename T_idx_>
struct Suffix_Array<T_idx_>::Subproblem_Ext_Mem
{
    Ext_Mem_Bucket<idx_t> SA_bucket;    // External-memory buckets for the SA elements.
    Ext_Mem_Bucket<idx_t> LCP_bucket;   // External-memory buckets for the LCP-array elements.
    Ext_Mem_Bucket<idx_t> sz_bucket;    // External-memory buckets for the sizes of the sorted sub-subarrays.


    Subproblem_Ext_Mem(const Suffix_Array<T_idx_>& SA, const std::size_t p_id):
          SA_bucket(SA.SA_bucket_file_path(p_id))
        , LCP_bucket(SA.LCP_bucket_file_path(p_id))
        , sz_bucket(SA.sz_bucket_file_path(p_id))
    {}

    Subproblem_Ext_Mem(Subproblem_Ext_Mem&&) = default;

    Subproblem_Ext_Mem(const Subproblem_Ext_Mem&) = delete;
    Subproblem_Ext_Mem& operator=(const Subproblem_Ext_Mem&) = delete;
    Subproblem_Ext_Mem&& operator=(Subproblem_Ext_Mem&&) = delete;

    void close()
    {
        SA_bucket.close(),
        LCP_bucket.close(),
        sz_bucket.close();
    }

    void remove()
    {
        SA_bucket.remove(),
        LCP_bucket.remove(),
        sz_bucket.remove();
    }
};


template <typename T_idx_>
inline T_idx_ Suffix_Array<T_idx_>::lcp(const char* const x, const char* const y, const idx_t min_len)
{
    idx_t l = 0;
    while(l < min_len && x[l] == y[l])
        l++;

    return l;
}


template <typename T_idx_>
template <std::size_t N>
inline T_idx_ Suffix_Array<T_idx_>::LCP(const char* const x, const char* const y, const idx_t min_len)
{
    idx_t lcp = 0;

    if constexpr(N == 1)
    {
        for(; lcp < min_len; ++lcp)
            if(x[lcp] != y[lcp])
                break;

        return lcp;
    }
    else
    {
        while((min_len - lcp) >= N * 32)
        {
            const auto l = LCP_unrolled<N>(x + lcp, y + lcp);
            lcp += l;
            if(l < N * 32)
                return lcp;
        }

        return lcp + LCP<N - 1>(x + lcp, y + lcp, min_len - lcp);
    }
}


template <typename T_idx_>
template <std::size_t N>
inline T_idx_ Suffix_Array<T_idx_>::LCP_unrolled(const char* const x, const char* const y)
{
    if constexpr(N == 0)
        return 0;
    else
    {
        const auto v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x));
        const auto v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y));
        const auto cmp = _mm256_cmpeq_epi8(v1, v2);
        const auto mask = static_cast<uint32_t>(_mm256_movemask_epi8(cmp));
        if(mask != 0xFFFFFFFF)
            return __builtin_ctz(~mask);

        return 32 + LCP_unrolled<N - 1>(x + 32, y + 32);
    }
}

}



#endif
