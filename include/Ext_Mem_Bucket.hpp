
#ifndef EXT_MEM_BUCKET_HPP
#define EXT_MEM_BUCKET_HPP



#include <cstddef>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <algorithm>
#include <cassert>


namespace CaPS_SA
{

// =============================================================================
// An external-memory-backed bucket for elements of type `T_`.
template <typename T_>
class Ext_Mem_Bucket
{
private:

    static constexpr std::size_t in_memory_bytes = 16lu * 1024; // 16KB.

    const std::string file_path;    // Path to the file storing the bucket.
    const std::size_t max_write_buf_bytes;  // Maximum size of the in-memory write-buffer in bytes.
    const std::size_t max_write_buf_elems;  // Maximum size of the in-memory write-buffer in elements.

    std::vector<T_> buf;    // In-memory buffer of the bucket-elements.
    std::size_t size_;  // Number of elements added to the bucket.

    std::ofstream file; // The bucket-file.


    // Flushes the in-memory buffer content to external memory.
    void flush();


public:

    // Constructs an external-memory bucket at path `file_path`. An optional in-
    // memory buffer size (in bytes) `buf_sz` for the bucket can be specified.
    // Specifying `append` indicates that the bucket exists and new elements are
    // to be appended to it.
    Ext_Mem_Bucket(const std::string& file_path, bool append = false, const std::size_t buf_sz = in_memory_bytes);

    // Move-constructs the bucket.
    Ext_Mem_Bucket(Ext_Mem_Bucket&& other) = default;

    // Invalidate the copy-constructor.
    Ext_Mem_Bucket(const Ext_Mem_Bucket& other) = delete;

    // Returns the size of the bucket.
    std::size_t size() const { return size_; }

    // Adds the element `elem` to the bucket.
    void add(const T_& elem);

    // Adds `sz` elements from `src` to the bucket.
    void add(const T_* src, std::size_t sz);

    // Dumps `sz` elements from `buf` directly into the external-memory. The
    // current in-memory elements of the bucket are bypassed.
    void dump(const T_* buf, std::size_t sz);

    // Closes the bucket. Elements should not be added anymore once this has
    // been invoked. This method is required only if the entirety of the bucket
    // needs to live in external-memory after the parent process finishes.
    void close();

    // Loads the bucket into the vector `v`.
    void load(std::vector<T_>& v) const;
};


template <typename T_>
inline Ext_Mem_Bucket<T_>::Ext_Mem_Bucket(const std::string& file_path, const bool append, const std::size_t buf_sz):
      file_path(file_path)
    , max_write_buf_bytes(buf_sz)
    , max_write_buf_elems(buf_sz / sizeof(T_))
    , size_(0)
{
    buf.reserve(max_write_buf_elems);

    if(append)
    {
        assert(!file_path.empty());

        std::error_code ec;
        const auto file_sz = std::filesystem::file_size(file_path, ec);
        if(ec)
        {
            std::cerr << "Error reading of external-memory bucket at " << file_path << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        assert(file_sz % sizeof(T_) == 0);
        size_ = file_sz / sizeof(T_);

        file.open(file_path, std::ios::app | std::ios::binary);
    }
    else if(!file_path.empty())
        file.open(file_path, std::ios::out | std::ios::binary);
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::add(const T_& elem)
{
    buf.push_back(elem);
    size_++;
    if(buf.size() >= max_write_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::add(const T_* const src, const std::size_t sz)
{
    while(buf.capacity() < buf.size() + sz)
        buf.reserve(buf.capacity() * 2);

    size_ += sz;
    std::copy(src, src + sz, buf.end());
    if(buf.size() >= max_write_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::dump(const T_* const buf, const std::size_t sz)
{
    size_ += sz;
    file.write(reinterpret_cast<const char*>(buf), sz * sizeof(T_));
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::flush()
{
    file.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
    buf.clear();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::close()
{
    if(!buf.empty())
        flush();

    file.close();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::load(std::vector<T_>& v) const
{
    std::error_code ec;
    const auto file_sz = std::filesystem::file_size(file_path, ec);

    assert(file_sz % sizeof(T_) == 0);
    assert(file_sz / sizeof(T_) + buf.size() == size_);

    v.resize(size_);

    std::ifstream input(file_path);
    input.read(reinterpret_cast<char*>(v.data()), file_sz);
    input.close();

    if(ec || !input)
    {
        std::cerr << "Error reading of external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::memcpy(reinterpret_cast<char*>(v.data()) + file_sz, reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
}

}



#endif
