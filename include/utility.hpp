
#ifndef CAPS_SA_UTILITY_HPP
#define CAPS_SA_UTILITY_HPP



#include <cstddef>
#include <utility>
#include <cstdlib>
#include <chrono>


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

}



#endif
