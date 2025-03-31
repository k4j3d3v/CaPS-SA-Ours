
#ifndef CAPS_SA_UTILITY_HPP
#define CAPS_SA_UTILITY_HPP



#include <utility>


namespace CaPS_SA
{

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
