
#ifndef EXT_MEM_BUCKET_HPP
#define EXT_MEM_BUCKET_HPP



#include <cstddef>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
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

    static constexpr std::size_t in_memory_bytes = 32lu * 1024; // 32KB.

    const std::string file_path;    // Path to the file storing the bucket.
    const std::size_t max_write_buf_bytes;  // Maximum size of the in-memory write-buffer in bytes.
    const std::size_t max_write_buf_elems;  // Maximum size of the in-memory write-buffer in elements.

    std::vector<T_> buf;    // In-memory buffer of the bucket-elements.
    std::size_t size_;  // Number of elements added to the bucket.

    const int fd;   // The bucket-file's descriptor.


    // Flushes the in-memory buffer content to external memory.
    void flush();

    // Checks the return code `code` of a file operation and reports pertinent
    // error(s) if any (also exits in that case). Otherwise returns `code`.
    int checked_file_op(const int code, const char* const op) const;

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
    , fd(checked_file_op(::open(file_path.c_str(), O_WRONLY | O_TRUNC | O_CREAT, S_IRWXU), "opening"))
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
inline void Ext_Mem_Bucket<T_>::flush()
{
    const auto r = checked_file_op(::write(fd, static_cast<const void*>(buf.data()), buf.size() * sizeof(T_)), " writing ");
    assert(static_cast<std::size_t>(r) == buf.size() * sizeof(T_)); (void)r;

    buf.clear();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::close()
{
    if(!buf.empty())
        flush();

    checked_file_op(::close(fd), " closing ");
}


template <typename T_>
inline std::size_t Ext_Mem_Bucket<T_>::load(T_* const dest) const
{
    const auto f = checked_file_op(::open(file_path.c_str(), O_RDONLY), " opening ");

    const auto to_read = size_ * sizeof(T_);
    const std::size_t max_read = 0x7ffff000;    // https://man7.org/linux/man-pages/man2/read.2.html#NOTES
    std::size_t read_bytes = 0;

    while(read_bytes < to_read)
    {
        const auto bytes = std::min(max_read, to_read - read_bytes);
        const auto r = checked_file_op(::read(f, reinterpret_cast<char*>(dest) + read_bytes, bytes), "reading");
        assert(static_cast<std::size_t>(r) == bytes); (void)r;

        read_bytes += r;
    }

    checked_file_op(::close(f), "closing");

    return size_;
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::rewrite(const T_* const src, const std::size_t sz)
{
    buf.clear();

    const auto f = checked_file_op(::open(file_path.c_str(), O_WRONLY | O_TRUNC), "opening");

    const auto to_write = sz * sizeof(T_);
    const std::size_t max_write = 0x7ffff000;    // https://man7.org/linux/man-pages/man2/write.2.html#NOTES
    std::size_t wrote_bytes = 0;

    while(wrote_bytes < to_write)
    {
        const auto bytes = std::min(max_write, to_write - wrote_bytes);
        const auto w = checked_file_op(::write(f, reinterpret_cast<const char*>(src) + wrote_bytes, bytes), "writing");
        assert(static_cast<std::size_t>(w) == bytes); (void)w;

        wrote_bytes += w;
    }

    checked_file_op(::close(f), "closing");

    size_ = sz;
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::remove()
{
    checked_file_op(::unlink(file_path.c_str()), "removing");
}


template <typename T_>
inline int Ext_Mem_Bucket<T_>::checked_file_op(const int code, const char* const op) const
{
    if(code == -1)
    {
        std::cerr << "Error " << op << " external-memory bucket at " << file_path << ". Aborting.\n";
        perror("Error");
        std::exit(EXIT_FAILURE);
    }

    return code;
}

}



#endif
