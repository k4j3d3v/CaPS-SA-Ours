
#ifndef CAPS_SA_UTILITY_HPP
#define CAPS_SA_UTILITY_HPP

#ifdef CAPSSA_QUIET
#define CAPS_SA_LOG(STATEMENT) ((void)0)
#else
#define CAPS_SA_LOG(STATEMENT) (STATEMENT)
#endif

#include <cstddef>
#include <utility>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <cassert>


namespace CaPS_SA
{

// TODO: merge the following two using nullptr `realloc`.

// Returns pointer to a memory-allocation for `size` elements of type `T_`.
template <typename T_>
static T_* allocate(std::size_t size) { return static_cast<T_*>(std::malloc(size * sizeof(T_))); }

// Deallocates the pointer `ptr`, allocated with `allocate`.
template <typename T_>
static void deallocate(T_* const ptr) { std::free(ptr); }

// Returns pointer to a memory-reallocation of pointer `ptr` for `size`
// elements of type `T_`.
template <typename T_>
static T_* reallocate(T_* const ptr, std::size_t size) { return static_cast<T_*>(std::realloc(ptr, size * sizeof(T_))); }

// Allocates the type-`T_` container `p` that currently has space for `cur_sz`
// elements geometrically with the growth factor `gf` such that it has enough
// space for at least `req_sz` elements, and returns the new size. If `keep_`
// is `true`, then the existing elements are kept.
template <typename T_, bool keep_ = false>
std::size_t reserve_geometric(T_*& p, const std::size_t curr_sz, const std::size_t req_sz, const double gf = 2.0)
{
    assert(gf > 1.0);

    if(curr_sz >= req_sz)
        return curr_sz;

    std::size_t new_sz = std::max(curr_sz, 1lu);
    while(new_sz < req_sz)
        new_sz *= gf;

    if constexpr(keep_)
        p = reallocate(p, new_sz);
    else
    {
        deallocate(p);
        p = allocate<T_>(new_sz);
    }

    return new_sz;
}

// Fields for profiling time.
typedef std::chrono::high_resolution_clock::time_point time_point_t;
constexpr static auto now = std::chrono::high_resolution_clock::now;
constexpr static auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


// Wrapper class for a data element of type `T_` to ensure that in a linear
// collection of `T_`'s, each element is aligned to a cache-line boundary.
template <typename T_>
class alignas(L1_CACHE_LINE_SIZE)
    Padded
{
private:

    T_ data_;


public:

    Padded()
    {}

    Padded(const T_& data):
      data_(data)
    {}

    Padded(T_&& data):
        data_(std::move(data))
    {}

    T_& unwrap() { return data_; }

    const T_& unwrap() const { return data_; }
};


// Wrapper class for a buffer of elements of type `T_`.
template <typename T_>
class Buffer
{
private:

    std::size_t cap_;   // Capacity of the buffer.
    T_* buf_;   // The raw buffer.


public:

    // Constructs a buffer with capacity `cap`.
    Buffer(const std::size_t cap):
          cap_(cap)
        , buf_(allocate<T_>(cap))
    {}

    ~Buffer() { deallocate(buf_); }

    Buffer(Buffer&& rhs) { *this = std::move(rhs); }

    Buffer& operator=(Buffer&& rhs)
    {
        cap_ = rhs.cap_;
        buf_ = rhs.buf_;
        rhs.cap_ = 0;
        rhs.buf_ = nullptr;

        return *this;
    }

    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;

    // Returns the memory region of the buffer.
    T_* data() { return buf_; }

    // Returns the memory region of the buffer.
    const T_* data() const { return buf_; }

    // Returns reference to the `idx`'th element of the buffer.
    T_& operator[](const std::size_t idx) { return buf_[idx]; }

    // Returns the `idx`'th element of the buffer.
    const T_& operator[](const std::size_t idx) const { return buf_[idx]; }

    // Returns the capacity of the buffer.
    auto capacity() const { return cap_; }

    // Ensures that the buffer have space for at least `new_cap` elements. No
    // guarantees are made for the existing elements.
    void reserve_uninit(const std::size_t new_cap) { cap_ = reserve_geometric(buf_, cap_, new_cap); }

    // Ensures that the buffer have space for at least `new_cap` elements.
    void reserve(const std::size_t new_cap) { cap_ = reserve_geometric<T_, true>(buf_, cap_, new_cap); }

    // Frees the buffer's memory.
    void free() { deallocate(buf_); buf_ = nullptr; cap_ = 0; }
};


// Returns the LCP length of `x` and `y`, with context-length `ctx`.
inline std::size_t lcp_unvectorized(const char* const x, const char* const y, const std::size_t ctx)
{
    std::size_t l = 0;
    while(l < ctx && x[l] == y[l])
        l++;

    return l;
}


template <typename T_seq_>
void read_input(const std::string& ip_path, std::vector<T_seq_>& text)
{
    std::error_code ec;
    const auto file_size = std::filesystem::file_size(ip_path, ec);

    if(ec)
    {
        CAPS_SA_LOG(std::cerr << ip_path << " : " << ec.message() << "\n");
        std::exit(EXIT_FAILURE);
    }

    assert(file_size % sizeof(T_seq_) == 0);

    text.resize(file_size / sizeof(T_seq_));
    std::ifstream input(ip_path);
    input.read(reinterpret_cast<char*>(text.data()), file_size);
    input.close();
}

}



#endif
