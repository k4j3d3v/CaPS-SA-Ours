
#ifndef EXT_MEM_BUCKET_HPP
#define EXT_MEM_BUCKET_HPP



#include <cstddef>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
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
    Ext_Mem_Bucket(const std::string& file_path, const std::size_t buf_sz = in_memory_bytes);

    // Move-constructs the bucket.
    Ext_Mem_Bucket(Ext_Mem_Bucket&& other) = default;

    // Invalidate the copy-constructor.
    Ext_Mem_Bucket(const Ext_Mem_Bucket&) = delete;

    // Invalidate the copy-assignment operator.
    Ext_Mem_Bucket& operator=(const Ext_Mem_Bucket&) = delete;

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

    // Loads the bucket into `dest`.
    std::size_t load(T_* dest) const;

    // Rewrites the bucket with `sz` elements from `src`.
    void rewrite(const T_* src, std::size_t sz);

    // Removes the bucket from disk.
    void remove();
};


template <typename T_>
inline Ext_Mem_Bucket<T_>::Ext_Mem_Bucket(const std::string& file_path, const std::size_t buf_sz):
      file_path(file_path)
    , max_write_buf_bytes(buf_sz)
    , max_write_buf_elems(buf_sz / sizeof(T_))
    , size_(0)
    , file(file_path, std::ios::out | std::ios::binary)
{
    assert(!file_path.empty());

    buf.reserve(max_write_buf_elems);
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
    size_ += sz;
    buf.insert(buf.end(), src, src + sz);
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
inline std::size_t Ext_Mem_Bucket<T_>::load(T_* const dest) const
{
    std::ifstream input(file_path);
    input.read(reinterpret_cast<char*>(dest), size_ * sizeof(T_));
    input.close();

    if(!input)
    {
        std::cerr << "Error reading of external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    assert(input.gcount() % sizeof(T_) == 0);
    return input.gcount() / sizeof(T_);
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::rewrite(const T_* const src, const std::size_t sz)
{
    buf.clear();

    if(!file.is_open())
        file.open(file_path, std::ios::out | std::ios::binary);

    file.clear();
    file.seekp(std::ios_base::beg);
    file.write(reinterpret_cast<const char*>(src), sz * sizeof(T_));

    size_ = sz;
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::remove()
{
    if(std::remove(file_path.c_str()))
    {
        std::cerr << "Error removing external-memory file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

}



#endif
