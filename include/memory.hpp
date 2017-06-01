#ifndef _MARRAY_MEMORY_HPP_
#define _MARRAY_MEMORY_HPP_

#include "utility.hpp"

namespace MArray
{

template <typename Type, typename Allocator>
class memory : private Allocator
{
    public:
        typedef std::allocator_traits<Allocator> alloc_traits;
        typedef alloc_traits::pointer pointer;
        typedef alloc_traits::const_pointer const_pointer;
        typedef alloc_traits::size_type size_type;

        memory() {}

        memory(const memory&) = delete;

        memory(memory&& other)
        : Allocator(std::move(other)), data_(other.data_), size_(other.size_)
        {
            other.data_ = nullptr;
            other.size_ = 0;
        }

        memory(const Allocator& alloc) : Allocator(alloc) {}

        ~memory()
        {
            destroy();
            deallocate();
        }

        memory& operator=(const memory&) = delete;

        memory& operator=(memory&& other)
        {
            Allocator::operator=(std::move(other));
            swap(other);
            return *this;
        }

        void allocate(size_type size)
        {
            data_ = alloc_traits::allocate(*this, size);
            size_ = (data_ ? size : 0);
        }

        void deallocate()
        {
            if (data_) alloc_traits::deallocate(*this, data_, size_);
            data_ = nullptr;
            size_ = 0;
        }

        void destroy()
        {
            for (size_type i = 0;i < size_;i++)
            {
                alloc_traits::destroy(*this, data_+i);
            }
        }

        operator const_pointer() const { return data_; }

        operator pointer() { return data_; }

        explicit operator bool() const { return data_; }

        void swap(memory& other)
        {
            using std::swap;
            swap(data_, other.data_);
            swap(size_, other.size_);
        }

        friend void swap(memory& a, memory& b)
        {
            a.swap(b);
        }

    protected:
        pointer data_ = nullptr;
        size_type size_ = 0;
};

}

#endif
