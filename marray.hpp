#ifndef _MARRAY_MARRAY_HPP_
#define _MARRAY_MARRAY_HPP_

#ifndef MARRAY_TEST
#define MARRAY_TEST(...)
#endif

#include <type_traits>
#include <array>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <utility>

#include "iterator.hpp"
#include "utility.hpp"

#define VECTOR_ALIGNMENT 16

#ifndef MARRAY_BASE_ALIGNMENT
#define MARRAY_BASE_ALIGNMENT 64
#endif

#ifndef MARRAY_STRIDE_ALIGNMENT
#define MARRAY_STRIDE_ALIGNMENT VECTOR_ALIGNMENT
#endif

namespace MArray
{
    template <typename T, size_t N> struct aligned_allocator
    {
        typedef T value_type;

        aligned_allocator() {}
        template <typename U, size_t M> aligned_allocator(const aligned_allocator<U, M>& other) {}

        T* allocate(size_t n)
        {
            void* ptr;
            int ret = posix_memalign(&ptr, N, n*sizeof(T));
            if (ret != 0) throw std::bad_alloc();
            return (T*)ptr;
        }

        void deallocate(T* ptr, size_t n)
        {
            free(ptr);
        }
    };

    template <typename T, size_t N, typename U, size_t M>
    bool operator==(const aligned_allocator<T, N>&, const aligned_allocator<U, M>&) { return true; }

    template <typename T, size_t N, typename U, size_t M>
    bool operator!=(const aligned_allocator<T, N>&, const aligned_allocator<U, M>&) { return false; }

    template <typename T, unsigned ndim> class const_marray;
    template <typename T, unsigned ndim> class marray;
    template <typename T, unsigned ndim, unsigned dim> class const_marray_ref;
    template <typename T, unsigned ndim, unsigned dim> class marray_ref;
    template <typename T, unsigned ndim, unsigned dim, unsigned newdim> class const_marray_slice;
    template <typename T, unsigned ndim, unsigned dim, unsigned newdim> class marray_slice;

    namespace slice
    {
        /*
         * The type all_t specifies a range [0,len_i) for an array
         * dimension i of length len_i (i.e. it selects all of the data along
         * that dimension).
         */
        struct all_t : public range_t<int>
        {
            constexpr all_t() : range_t<int>(0, 0) {}
        };

        constexpr all_t all;
    }

    /*
     * Represents a part of an array, where the first dim-1 out of ndim
     * dimensions have been indexed into. This type may be implicity converted
     * to an array view or further indexed. This particular reference type
     * is explicitly const, and may only be used to read the original array.
     */
    template <typename T, unsigned ndim, unsigned dim>
    class const_marray_ref
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            typedef typename const_marray<T, ndim>::size_type size_type;
            typedef typename const_marray<T, ndim>::idx_type idx_type;
            typedef typename const_marray<T, ndim>::value_type value_type;
            typedef typename const_marray<T, ndim>::pointer pointer;
            typedef typename const_marray<T, ndim>::const_pointer const_pointer;
            typedef typename const_marray<T, ndim>::reference reference;
            typedef typename const_marray<T, ndim>::const_reference const_reference;

        protected:
            const const_marray<T, ndim>& array_;
            size_type idx;

            const_marray_ref(const const_marray_ref& other) = default;

            const_marray_ref(const const_marray<T, ndim>& array_, size_type idx, idx_type i)
            : array_(array_), idx(idx+i*array_.stride_[dim-2]) {}

        public:
            const_marray_ref& operator=(const const_marray_ref& other) = delete;

            const_marray_ref<T, ndim, dim+1> operator[](idx_type i) const
            {
                assert(i >= 0 && i < array_.len_[dim-1]);
                return const_marray_ref<T, ndim, dim+1>(array_, idx, i);
            }

            template <typename I>
            const_marray_slice<T, ndim, dim+1, 1> operator[](const range_t<I>& x) const
            {
                return const_marray_slice<T, ndim, dim+1, 1>(array_, idx, {}, {}, x);
            }

            const_marray_slice<T, ndim, dim+1, 1> operator[](const slice::all_t& x) const
            {
                return const_marray_slice<T, ndim, dim+1, 1>(array_, idx, {}, {}, range(idx_type(), array_.len_[dim-1]));
            }

            const_pointer data() const
            {
                return array_.data_+idx;
            }
    };


    /*
     * Identical to const_marray_ref, except that the reference is no longer
     * const; it may be used to alter data owned by the original array. This
     * reference may only be constructed from non-const arrays or references.
     */
    template <typename T, unsigned ndim, unsigned dim>
    class marray_ref : public const_marray_ref<T, ndim, dim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            using typename const_marray_ref<T, ndim, dim>::size_type;
            using typename const_marray_ref<T, ndim, dim>::idx_type;
            using typename const_marray_ref<T, ndim, dim>::value_type;
            using typename const_marray_ref<T, ndim, dim>::pointer;
            using typename const_marray_ref<T, ndim, dim>::const_pointer;
            using typename const_marray_ref<T, ndim, dim>::reference;
            using typename const_marray_ref<T, ndim, dim>::const_reference;

        protected:
            using const_marray_ref<T, ndim, dim>::array_;
            using const_marray_ref<T, ndim, dim>::idx;

            marray_ref(const marray_ref& other) = default;

            marray_ref(const const_marray<T, ndim>& array_, size_type idx, idx_type i)
            : const_marray_ref<T, ndim, dim>(array_, idx, i) {}

        public:
            marray_ref& operator=(const marray_ref& other)
            {
                copy(other, *this);
                return *this;
            }

            template <unsigned ndim_>
            marray_ref& operator=(const const_marray_ref<T, ndim_, ndim_-ndim+dim>& other)
            {
                copy(other, *this);
                return *this;
            }

            template <unsigned ndim_, unsigned newdim_>
            marray_ref& operator=(const const_marray_slice<T, ndim_, ndim_-ndim+dim+newdim_, newdim_>& other)
            {
                copy(other, *this);
                return *this;
            }

            marray_ref<T, ndim, dim+1> operator[](idx_type i)
            {
                assert(i >= 0 && i < array_.len_[dim-1]);
                return marray_ref<T, ndim, dim+1>(array_, idx, i);
            }

            template <typename I>
            marray_slice<T, ndim, dim+1, 1> operator[](const range_t<I>& x)
            {
                return marray_slice<T, ndim, dim+1, 1>(array_, idx, {}, {}, x);
            }

            marray_slice<T, ndim, dim+1, 1> operator[](const slice::all_t& x)
            {
                return marray_slice<T, ndim, dim+1, 1>(array_, idx, {}, {}, range(idx_type(), array_.len_[ndim-1]));
            }

            pointer data()
            {
                return array_.data_+idx;
            }
    };

    /*
     * Represents a part of an array, where all but the last dimension have
     * been indexed into. May be converted into an array view or indexed one
     * last time to produce a reference to a particular data element. This
     * version is explicitly const.
     */
    template <typename T, unsigned ndim>
    class const_marray_ref<T, ndim, ndim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            typedef typename marray<T, ndim>::size_type size_type;
            typedef typename marray<T, ndim>::idx_type idx_type;
            typedef typename marray<T, ndim>::value_type value_type;
            typedef typename marray<T, ndim>::pointer pointer;
            typedef typename marray<T, ndim>::const_pointer const_pointer;
            typedef typename marray<T, ndim>::reference reference;
            typedef typename marray<T, ndim>::const_reference const_reference;

        protected:
            const const_marray<T, ndim>& array_;
            size_type idx;

            const_marray_ref(const const_marray_ref& other) = default;

            const_marray_ref(const const_marray<T, ndim>& array_, size_type idx, idx_type i)
            : array_(array_), idx(idx+i*array_.stride_[ndim-2]) {}

        public:
            const_marray_ref& operator=(const const_marray_ref& other) = delete;

            const_reference operator[](idx_type i) const
            {
                assert(i >= 0 && i < array_.len_[ndim-1]);
                return array_.data_[idx+i*array_.stride_[ndim-1]];
            }

            template <typename I>
            const_marray<T, 1> operator[](const range_t<I>& x) const
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array_.len_[ndim-1]);
                const_marray<T, 1> ret(*this);
                ret.size_ = x.size();
                ret.data_ += ret.stride_*x.front();
                return ret;
            }

            const_marray<T, 1> operator[](const slice::all_t& x) const
            {
                return const_marray<T, 1>(*this);
            }

            const_pointer data() const
            {
                return array_.data_+idx;
            }
    };


    /*
     * Identical to const_marray_ref<T,ndim,ndim> except that it is non-const.
     */
    template <typename T, unsigned ndim>
    class marray_ref<T, ndim, ndim> : public const_marray_ref<T, ndim, ndim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            using typename const_marray_ref<T, ndim, ndim>::size_type;
            using typename const_marray_ref<T, ndim, ndim>::idx_type;
            using typename const_marray_ref<T, ndim, ndim>::value_type;
            using typename const_marray_ref<T, ndim, ndim>::pointer;
            using typename const_marray_ref<T, ndim, ndim>::const_pointer;
            using typename const_marray_ref<T, ndim, ndim>::reference;
            using typename const_marray_ref<T, ndim, ndim>::const_reference;

        protected:
            using const_marray_ref<T, ndim, ndim>::array_;
            using const_marray_ref<T, ndim, ndim>::idx;

            marray_ref(const marray_ref& other) = default;

            marray_ref(const const_marray<T, ndim>& array_, size_type idx, idx_type i)
            : const_marray_ref<T, ndim, ndim>(array_, idx, i) {}

        public:
            marray_ref& operator=(const marray_ref& other)
            {
                copy(other, *this);
                return *this;
            }

            template <unsigned ndim_>
            marray_ref& operator=(const const_marray_ref<T, ndim_, ndim_>& other)
            {
                copy(other, *this);
                return *this;
            }

            reference operator[](idx_type i)
            {
                assert(i >= 0 && i < array_.len_[ndim-1]);
                return const_cast<reference>(array_.data_[idx+i*array_.stride_[ndim-1]]);
            }

            template <typename I>
            marray<T, 1> operator[](const range_t<I>& x)
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array_.len_[ndim-1]);
                marray<T, 1> ret(*this);
                ret.size_ = x.size();
                ret.data_ += ret.stride_*x.front();
                return ret;
            }

            marray<T, 1> operator[](const slice::all_t& x)
            {
                return marray<T, 1>(*this);
            }

            pointer data()
            {
                return array_.data_+idx;
            }
    };

    /*
     * Represents a part of an array, where the first dim-1 out of ndim
     * dimensions have either been indexed into (i.e. a single value
     * specified for that index) or sliced (i.e. a range of values specified).
     * The parameter newdim specifies how many indices were sliced. The
     * reference may be converted into an array view (of dimension
     * ndim-dim+1+newdim) or further indexed, but may not be used to modify
     * data.
     */
    template <typename T, unsigned ndim, unsigned dim, unsigned newdim>
    class const_marray_slice
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            typedef typename marray<T, ndim>::size_type size_type;
            typedef typename marray<T, ndim>::idx_type idx_type;
            typedef typename marray<T, ndim>::value_type value_type;
            typedef typename marray<T, ndim>::pointer pointer;
            typedef typename marray<T, ndim>::const_pointer const_pointer;
            typedef typename marray<T, ndim>::reference reference;
            typedef typename marray<T, ndim>::const_reference const_reference;

        protected:
            const const_marray<T, ndim>& array_;
            size_type idx;
            std::array<unsigned, newdim> dims;
            std::array<idx_type, newdim> lens;

            const_marray_slice(const const_marray_slice& other) = default;

            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               const std::array<unsigned,newdim>& dims,
                               const std::array<idx_type,newdim>& lens, idx_type i)
            : array_(array_), idx(idx+i*array_.stride[dim-2]), dims(dims), lens(lens) {}

            template <typename I>
            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               const std::array<unsigned,newdim-1>& dims,
                               const std::array<idx_type,newdim-1>& lens, const range_t<I>& range_)
            : array_(array_), idx(idx+array_.stride_[dim-2]*range_.front())
            {
                std::copy(dims.begin(), dims.end(), this->dims.begin());
                this->dims.back() = dim-2;
                std::copy(lens.begin(), lens.end(), this->lens.begin());
                this->lens.back() = range_.size();
            }

        public:
            const_marray_slice& operator=(const const_marray_slice& other) = delete;

            const_marray_slice<T, ndim, dim+1, newdim> operator[](idx_type i) const
            {
                assert(i >= 0 && i < array_.len[dim-1]);
                return const_marray_slice<T, ndim, dim+1, newdim>(array_, idx, dims, lens, i);
            }

            template <typename I>
            const_marray_slice<T, ndim, dim+1, newdim+1> operator[](const range_t<I>& x) const
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array_.len[dim-1]);
                return const_marray_slice<T, ndim, dim+1, newdim+1>(array_, idx, dims, lens, x);
            }

            const_marray_slice<T, ndim, dim+1, newdim+1> operator[](const slice::all_t& x) const
            {
                return const_marray_slice<T, ndim, dim+1, newdim+1>(array_, idx, dims, lens, range(idx_type(), array_.len_[dim-1]));
            }

            const_pointer data() const
            {
                return array_.data_+idx;
            }
    };

    /*
     * Identical to const_marray_slice, except that this version is non-const
     * and may be used to modify data.
     */
    template <typename T, unsigned ndim, unsigned dim, unsigned newdim>
    class marray_slice : public const_marray_slice<T, ndim, dim, newdim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            using typename const_marray_slice<T, ndim, dim, newdim>::size_type;
            using typename const_marray_slice<T, ndim, dim, newdim>::idx_type;
            using typename const_marray_slice<T, ndim, dim, newdim>::value_type;
            using typename const_marray_slice<T, ndim, dim, newdim>::pointer;
            using typename const_marray_slice<T, ndim, dim, newdim>::const_pointer;
            using typename const_marray_slice<T, ndim, dim, newdim>::reference;
            using typename const_marray_slice<T, ndim, dim, newdim>::const_reference;

        protected:
            using const_marray_slice<T, ndim, dim, newdim>::array_;
            using const_marray_slice<T, ndim, dim, newdim>::idx;
            using const_marray_slice<T, ndim, dim, newdim>::dims;
            using const_marray_slice<T, ndim, dim, newdim>::lens;

            marray_slice(const marray_slice& other) = default;

            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                         const std::array<unsigned,newdim>& dims,
                         const std::array<idx_type,newdim>& lens, idx_type i)
            : const_marray_slice<T, ndim, dim, newdim>(array_, idx, dims, lens, i) {}

            template <typename I>
            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                         const std::array<unsigned,newdim-1>& dims,
                         const std::array<idx_type,newdim-1>& lens, const range_t<I>& range_)
            : const_marray_slice<T, ndim, dim, newdim>(array_, idx, dims, lens, range_) {}

        public:
            marray_slice& operator=(const marray_slice& other)
            {
                copy(other, *this);
                return *this;
            }

            template <unsigned ndim_>
            marray_slice& operator=(const const_marray_ref<T, ndim_, ndim_-ndim+dim-newdim>& other)
            {
                copy(other, *this);
                return *this;
            }

            template <unsigned ndim_, unsigned newdim_>
            marray_slice& operator=(const const_marray_slice<T, ndim_, ndim_-ndim+dim-newdim+newdim_, newdim_>& other)
            {
                copy(other, *this);
                return *this;
            }

            marray_slice<T, ndim, dim+1, newdim> operator[](idx_type i)
            {
                assert(i >= 0 && i < array_.len_[dim-1]);
                return marray_slice<T, ndim, dim+1, newdim>(array_, idx, dims, lens, i);
            }

            template <typename I>
            marray_slice<T, ndim, dim+1, newdim+1> operator[](const range_t<I>& x)
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array_.len_[dim-1]);
                return marray_slice<T, ndim, dim+1, newdim+1>(array_, idx, dims, lens, x);
            }

            marray_slice<T, ndim, dim+1, newdim+1> operator[](const slice::all_t& x)
            {
                return marray_slice<T, ndim, dim+1, newdim+1>(array_, idx, dims, lens, range(idx_type(), array_.len_[ndim-1]));
            }

            pointer data()
            {
                return array_.data_+idx;
            }
    };

    /*
     * Represents a part of an array, where all but the last of ndim dimensions
     * have either been indexed into (i.e. a single value specified for that
     * index) or sliced (i.e. a range of values specified). The parameter newdim
     * specifies how many indices were sliced. The reference may be converted
     * into an array view (of dimension newdim+1) or indexed one last
     * time to create a view of dimension newdim (if a single value is specified)
     * or newdim+1 (if a range is specified). This version is explicitly const.
     */
    template <typename T, unsigned ndim, unsigned newdim>
    class const_marray_slice<T, ndim, ndim, newdim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            typedef typename marray<T, ndim>::size_type size_type;
            typedef typename marray<T, ndim>::idx_type idx_type;
            typedef typename marray<T, ndim>::value_type value_type;
            typedef typename marray<T, ndim>::pointer pointer;
            typedef typename marray<T, ndim>::const_pointer const_pointer;
            typedef typename marray<T, ndim>::reference reference;
            typedef typename marray<T, ndim>::const_reference const_reference;

        protected:
            const const_marray<T, ndim>& array_;
            size_type idx;
            std::array<unsigned, newdim> dims;
            std::array<idx_type, newdim> lens;

            const_marray_slice(const const_marray_slice& other) = default;

            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               const std::array<unsigned,newdim>& dims,
                               const std::array<idx_type,newdim>& lens, idx_type i)
            : array_(array_), idx(idx+i*array_.stride_[ndim-2]), dims(dims), lens(lens) {}

            template <typename I>
            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               const std::array<unsigned,newdim-1>& dims,
                               const std::array<idx_type,newdim-1>& lens, const range_t<I>& x)
            : array_(array_), idx(idx+x.front()*array_.stride_[ndim-2])
            {
                std::copy(dims.begin(), dims.end(), this->dims.begin());
                this->dims.back() = ndim-2;
                std::copy(lens.begin(), lens.end(), this->lens.begin());
                this->lens.back() = x.size();
            }

        public:
            const_marray_slice& operator=(const const_marray_slice& other) = delete;

            const_marray<T, newdim> operator[](idx_type i) const
            {
                assert(i >= 0 && i < array_.len_[ndim-1]);
                const_marray<T, newdim> ret;
                ret.data_ = array_.data_+idx+i*array_.stride_[ndim-1];
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = lens[dim];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                return ret;
            }

            template <typename I>
            const_marray<T, newdim+1> operator[](const range_t<I>& x) const
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array_.len_[ndim-1]);
                const_marray<T, newdim+1> ret;
                ret.data_ = array_.data_+idx+array_.stride_[ndim-1]*x.front();
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = lens[dim];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                ret.len_[newdim] = x.size();
                ret.stride_[newdim] = array_.stride_[ndim-1];
                ret.size_ *= ret.len_[newdim];
                return ret;
            }

            const_marray<T, newdim+1> operator[](const slice::all_t& x) const
            {
                return const_marray<T, newdim+1>(*this);
            }

            const_pointer data() const
            {
                return array_.data_+idx;
            }
    };


    /*
     * Identical to const_marray_slice except that the data may be modified.
     */
    template <typename T, unsigned ndim, unsigned newdim>
    class marray_slice<T, ndim, ndim, newdim> : public const_marray_slice<T, ndim, ndim, newdim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            using typename const_marray_slice<T, ndim, ndim, newdim>::size_type;
            using typename const_marray_slice<T, ndim, ndim, newdim>::idx_type;
            using typename const_marray_slice<T, ndim, ndim, newdim>::value_type;
            using typename const_marray_slice<T, ndim, ndim, newdim>::pointer;
            using typename const_marray_slice<T, ndim, ndim, newdim>::const_pointer;
            using typename const_marray_slice<T, ndim, ndim, newdim>::reference;
            using typename const_marray_slice<T, ndim, ndim, newdim>::const_reference;

        protected:
            using const_marray_slice<T, ndim, ndim, newdim>::array_;
            using const_marray_slice<T, ndim, ndim, newdim>::idx;
            using const_marray_slice<T, ndim, ndim, newdim>::dims;
            using const_marray_slice<T, ndim, ndim, newdim>::lens;

            marray_slice(const marray_slice& other) = default;

            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                         const std::array<unsigned,newdim>& dims,
                         const std::array<idx_type,newdim>& lens, idx_type i)
            : const_marray_slice<T, ndim, ndim, newdim>(array_, idx, dims, lens, i) {}

            template <typename I>
            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                         const std::array<unsigned,newdim-1>& dims,
                         const std::array<idx_type,newdim-1>& lens, const range_t<I>& x)
            : const_marray_slice<T, ndim, ndim, newdim>(array_, idx, dims, lens, x) {}

        public:
            marray_slice& operator=(const marray_slice& other)
            {
                copy(other, *this);
                return *this;
            }

            template <unsigned ndim_>
            marray_slice& operator=(const const_marray_ref<T, ndim_, ndim_-newdim>& other)
            {
                copy(other, *this);
                return *this;
            }

            template <unsigned ndim_, unsigned newdim_>
            marray_slice& operator=(const const_marray_slice<T, ndim_, ndim_-newdim+newdim_, newdim_>& other)
            {
                copy(other, *this);
                return *this;
            }

            marray<T, newdim> operator[](idx_type i)
            {
                assert(i >= 0 && i < array_.len_[ndim-1]);
                marray<T, newdim> ret;
                ret.data_ = array_.data_+idx+i*array_.stride_[ndim-1];
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = lens[dim];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                return ret;
            }

            template <typename I>
            marray<T, newdim+1> operator[](const range_t<I>& x)
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array_.len_[ndim-1]);
                marray<T, newdim+1> ret;
                ret.data_ = array_.data_+idx+array_.stride_[ndim-1]*x.front();
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = lens[dim];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                ret.len_[newdim] = x.size();
                ret.stride_[newdim] = array_.stride_[ndim-1];
                ret.size_ *= ret.len_[newdim];
                return ret;
            }

            marray<T, newdim+1> operator[](const slice::all_t& x)
            {
                return marray<T, newdim+1>(*this);
            }

            pointer data()
            {
                return array_.data_+idx;
            }
    };

    /*
     * The special value construct_view is used to create a view of an array,
     * even if the default copy constructor would have created a copy of the
     * data.
     */
    struct construct_view_t {};
    constexpr construct_view_t construct_view = construct_view_t();

    /*
     * The special value construct_copy is used to create a copy of an array
     * or view, even if the default copy constructor would have created a view.
     */
    struct construct_copy_t {};
    constexpr construct_copy_t construct_copy = construct_copy_t();

    /*
     * The special value uninitialized is used to construct an array which
     * does not default- or value-initialize its elements (useful for avoiding
     * redundant memory operations for scalar types).
     */
    struct uninitialized_t {};
    constexpr uninitialized_t uninitialized = uninitialized_t();

    /*
     * Specifies the layout of the array data.
     */
    enum Layout : int {COLUMN_MAJOR, ROW_MAJOR, DEFAULT=ROW_MAJOR};

    namespace detail
    {
        inline size_t align(size_t n, size_t alignment)
        {
            return ((n+alignment-1)/alignment)*alignment;
        }

        /*
         * These helper structs determine the type of an marray with
         * datatype T, number of dimensions N, and const iff isConst == true.
         * Additionally, the N=0 case returns the type of a const or non-const
         * reference to the underlying datatype.
         */
        template <bool isConst, typename T, unsigned N>
        struct marray_type { typedef marray<T, N> type; };
        
        template <typename T, unsigned N>
        struct marray_type<true, T, N> { typedef const_marray<T, N> type; };

        template <typename T>
        struct marray_type<false, T, 0> { typedef typename marray<T, 1>::reference type; };

        template <typename T>
        struct marray_type<true, T, 0> { typedef typename marray<T, 1>::const_reference type; };

        /*
         * These helper structs determine the number of sliced dimensions based
         * on a set of indices which may be either an integral type (no slicing)
         * or a range (slicing).
         */
        template <typename Index>
        struct is_slice { static constexpr unsigned N = 0; };

        template <typename I>
        struct is_slice<range_t<I>> { static constexpr unsigned N = 1; };

        template <typename Index, typename... Indices>
        struct num_slices
        {
                static constexpr unsigned N = num_slices<Indices...>::N + is_slice<Index>::N;
        };

        template <typename Index>
        struct num_slices<Index>
        {
                static constexpr unsigned N = is_slice<Index>::N;
        };

        /*
         * This helper struct determines the return type after indexing an
         * array with the given index types. The return type may be either
         * a reference to the underlying datatype (no slicing), or a lower-
         * dimensional view of the array if one or more arguments are ranges.
         */
        template <bool isConst, typename T, typename... Indices>
        struct return_type
        {
            typedef typename marray_type<isConst, T, num_slices<Indices...>::N>::type type;
        };

        /*
         * These helper structs apply the given set of indices to an array
         * to generate either a reference to a specific element or a
         * lower-dimensional view of the array.
         */
        template <typename RT, typename Index, typename... Indices>
        struct get_slice
        {
            template <typename ArrayOrSlice>
            RT operator()(ArrayOrSlice&& x, const Index& idx0, const Indices&... idx)
            {
                return get_slice<RT, Indices...>()(x[idx0], idx...);
            }
        };

        template <typename RT, typename Index>
        struct get_slice<RT, Index>
        {
            template <typename ArrayOrSlice>
            RT operator()(ArrayOrSlice&& x, const Index& idx0)
            {
                return x[idx0];
            }
        };

        template <unsigned N, unsigned I, template <typename...> class Condition, typename... Args>
        struct are_integral_ : std::false_type {};

        template <unsigned N, unsigned I, template <typename...> class Condition, typename Arg, typename... Args>
        struct are_integral_<N, I, Condition, Arg, Args...>
        : std::conditional<std::is_integral<Arg>::value,
                           are_integral_<N, I+1, Condition, Args...>,
                           std::false_type>::type {};

        template <unsigned N, template <typename...> class Condition, typename Arg, typename... Args>
        struct are_integral_<N, N, Condition, Arg, Args...>
        : std::conditional<std::is_integral<Arg>::value,
                           Condition<Args...>,
                           std::false_type>::type {};

        template <unsigned N, template <typename...> class Condition, typename... Args>
        struct are_integral : are_integral_<N, 1, Condition, Args...> {};

        template <typename... Args>
        struct are_empty : std::integral_constant<bool, (sizeof...(Args) == 0)> {};

        /*
         * Helper classes to set the len_ and stride_ arrays of an marray.
         */
        template <typename T, unsigned N, unsigned I, typename Tail, typename Arg, typename... Args>
        struct set_len_
        {
            set_len_(const_marray<T, N>& array, std::array<size_t, N>& len, const Arg& arg, const Args&... args)
            {
                len[I-1] = arg;
                set_len_<T, N, I+1, Tail, Args...> tmp(array, len, args...);
            }
        };

        template <typename T, unsigned N, typename Tail, typename Arg, typename... Args>
        struct set_len_<T, N, N, Tail, Arg, Args...>
        {
            set_len_(const_marray<T, N>& array, std::array<size_t, N>& len, const Arg& arg, const Args&... args)
            {
                len[N-1] = arg;
                Tail tmp(array, len, args...);
            }
        };

        template <typename T, unsigned N, typename Tail, typename... Args>
        using set_len = set_len_<T, N, 1, Tail, Args...>;

        /*
         * Helper class to determine if the arguments Args... are of one of the
         * forms which are acceptable for constructing an instance of
         * [const_]marray<T, ndim> or passing to marray<T, ndim>::reset.
         */
        template <typename T, unsigned ndim, typename... Args>
        struct are_reset_args
        {
            template <typename... Args_>
            struct are_direct_args : std::true_type {};

            template <typename... Args_>
            struct are_ptr_args : std::false_type {};

            template <typename... Args_>
            struct are_remaining_args
            : std::integral_constant<bool, are_direct_args<Args_...>::value ||
                                           are_ptr_args<Args_...>::value> {};

            constexpr static bool value = are_integral<ndim, are_remaining_args, typename std::decay<Args>::type...>::value;

            struct reset_
            {
                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, Layout layout = DEFAULT)
                {
                    array.reset(len, layout);
                }

                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, const T& val, Layout layout = DEFAULT)
                {
                    array.reset(len, val, layout);
                }

                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, const uninitialized_t& u, Layout layout = DEFAULT)
                {
                    array.reset(len, u, layout);
                }

                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, unsigned align_at, Layout layout = DEFAULT)
                {
                    array.reset(len, align_at, layout);
                }

                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, const T& val, unsigned align_at, Layout layout = DEFAULT)
                {
                    array.reset(len, val, align_at, layout);
                }

                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
                {
                    array.reset(len, u, align_at, layout);
                }

                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, const T* ptr, Layout layout = DEFAULT)
                {
                    array.reset(len, const_cast<T*>(ptr), layout);
                }

                template <typename... Args_>
                reset_(const_marray<T, ndim>& array, std::array<size_t, ndim>& len, const T* ptr, const Args_&... args)
                {
                    array.reset(len, const_cast<T*>(ptr), make_array(args...));
                }
            };

            are_reset_args(const_marray<T, ndim>& array, const Args&... args)
            {
                std::array<size_t, ndim> len;
                set_len<T, ndim, reset_, Args...> tmp(array, len, args...);
            }
        };

        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg>
        struct are_reset_args<T, ndim, Args...>::are_direct_args<Arg>
        : std::integral_constant<bool, std::is_convertible<Arg,T>::value ||
                                       std::is_same<Arg,Layout>::value ||
                                       std::is_same<Arg,uninitialized_t>::value ||
                                       std::is_integral<Arg>::value> {};

        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg1, typename Arg2>
        struct are_reset_args<T, ndim, Args...>::are_direct_args<Arg1, Arg2>
        : std::integral_constant<bool, ((std::is_convertible<Arg1,T>::value ||
                                         std::is_same<Arg1,uninitialized_t>::value) &&
                                        (std::is_same<Arg2,Layout>::value ||
                                         std::is_integral<Arg2>::value)) ||
                                       ( std::is_integral<Arg1>::value &&
                                         std::is_same<Arg2,Layout>::value)> {};

        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg1, typename Arg2, typename Arg3>
        struct are_reset_args<T, ndim, Args...>::are_direct_args<Arg1, Arg2, Arg3>
        : std::integral_constant<bool, (std::is_convertible<Arg1,T>::value ||
                                        std::is_same<Arg1,uninitialized_t>::value) &&
                                        std::is_integral<Arg2>::value &&
                                        std::is_same<Arg3,Layout>::value> {};

        /*
         * This specialization is necessary to work around (different!) compiler
         * bugs in GCC and Clang.
         */
        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename... Args_>
        struct are_reset_args<T, ndim, Args...>::are_direct_args<Arg1, Arg2, Arg3, Arg4, Args_...> : std::false_type {};

        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg>
        struct are_reset_args<T, ndim, Args...>::are_ptr_args<Arg> : std::is_convertible<Arg,const T*> {};

        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg>
        struct are_reset_args<T, ndim, Args...>::are_ptr_args<Arg, Layout> : std::is_convertible<Arg,const T*> {};

        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg, typename... Args_>
        struct are_reset_args<T, ndim, Args...>::are_ptr_args<Arg, Args_...>
        : std::integral_constant<bool, std::is_convertible<Arg,const T*>::value &&
                                       are_integral<ndim, are_empty, Args_...>::value> {};

        /*
         * Helper class to determine if the arguments Args... are of one of the
         * forms which are acceptable for passing to marray<T, ndim>::resize.
         */
        template <typename T, unsigned ndim, typename... Args>
        struct are_resize_args
        {
            template <typename... Args_>
            struct are_direct_args : std::false_type {};

            template <typename... Args_>
            struct are_remaining_args
            : std::integral_constant<bool, are_empty<Args_...>::value ||
                                           are_direct_args<Args_...>::value> {};

            constexpr static bool value = are_integral<ndim, are_remaining_args, typename std::decay<Args>::type...>::value;

            struct resize_
            {
                resize_(const_marray<T, ndim>& array_, std::array<size_t, ndim>& len, const T& val = T())
                {
                    /*
                     * I'm too lazy to make set_len_ pass through a non-const reference.
                     */
                    marray<T, ndim>& array = static_cast<marray<T, ndim>&>(array_);
                    array.resize(len, val);
                }
            };

            template <typename... Args_>
            are_resize_args(marray<T, ndim>& array, const Args_&... args)
            {
                std::array<size_t, ndim> len;
                set_len<T, ndim, resize_, Args_...> tmp(array, len, args...);
            }
        };

        template <typename T, unsigned ndim, typename... Args>
        template <typename Arg>
        struct are_resize_args<T, ndim, Args...>::are_direct_args<Arg>
        : std::integral_constant<bool, std::is_convertible<Arg,T>::value> {};
    }

    /*
     * This class represents a multidimensional array of elements of type T,
     * with ndim > 0 dimensions. An array may either own its own data (in
     * which case it is responsible for allocation and deallocation), or
     * represent a "view" of all or part of another array. It is the user's
     * responsibility to ensure that views are not referenced once the original
     * array has been destructed or reset with the resize() or reset() methods,
     * since no internal reference counting is performed.
     *
     * Arrays are constructed using one of the following forms:
     *
     * 1) const_marray<T, ndim>()
     *
     *  Constructs an array with all dimensions of zero length.
     *
     * 2) const_marray<T, ndim>(const const_marray<T, ndim>& other)
     *
     *  Creates a copy of the given array. If the array is a view, then so is
     *  the copy.
     *
     * 3) const_marray<T, ndim>(const_marray<T, ndim>&& other)
     *
     *  Move-constructs a copy of the given array. If the array is a view, then
     *  this behaves identically to the copy constructor. Otherwise, the given
     *  array is reinitialized to a state as if it were constructed using the
     *  empty constructor.
     *
     * 4) const_marray<T, ndim>(const const_marray<T, ndim>& other,
     *                          const construct_view_t& cv)
     *
     *  Constructs a view of the given array, even if it is not a view itself.
     *
     * 5) const_marray<T, ndim>(const const_marray<T, ndim>& other,
     *                          const construct_copy_t& cp)
     *
     *  Constructs a new (freshly allocated) copy of the given array, even if
     *  it is a view.
     *
     * 6.1) template <typename U>
     *      explicit const_marray<T, ndim>(const std::array<U, ndim>& len,
     *                                     Layout layout = DEFAULT)
     *
     * 6.2) template <typename U>
     *      const_marray<T, ndim>(const std::array<U, ndim>& len,
     *                            const T& val, Layout layout = DEFAULT)
     *
     * 6.3) template <typename U>
     *      const_marray<T, ndim>(const std::array<U, ndim>& len,
     *                            const uninitialized_t& u, Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout. In the first
     *  form the elements are default-initialized. In the second form they are
     *  copy-initialized from val. In the third form they are not initialized.
     *
     * 7.4) template <typename U>
     *      explicit const_marray<T, ndim>(const std::array<U, ndim>& len,
     *                                     unsigned align_at, Layout layout = DEFAULT)
     *
     * 7.5) template <typename U>
     *      const_marray<T, ndim>(const std::array<U, ndim>& len,
     *                            const T& val, unsigned align_at, Layout layout = DEFAULT)
     *
     * 7.6) template <typename U>
     *      const_marray<T, ndim>(const std::array<U, ndim>& len, const uninitialized_t& u,
     *                            unsigned align_at, Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout. In the first
     *  form the elements are default-initialized. In the second form they are
     *  copy-initialized from val. In the third form they are not initialized.
     *  The strides of the indices are adjusted such that stride (align_at-1)
     *  for row-major or stride (align_at) for column-major storage is a
     *  multiple of the processor's vector length. align_at = 0 specifies no
     *  additional alignment. In all cases (even in constructors without an
     *  align_at parameter), the base pointer of the array is aligned to a
     *  multiple of the cache line size.
     *
     * 8.1) template <typename I1, ..., typename IN>
     *      explicit const_marray<T, ndim>(I1 len1, ..., IN lenN, Layout layout = DEFAULT)
     *
     * 8.2) template <typename I1, ..., typename IN>
     *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const T& val, Layout layout = DEFAULT)
     *
     * 8.3) template <typename I1, ..., typename IN>
     *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const unsigned_t& u,
     *                            Layout layout = DEFAULT)
     *
     * 9.1) template <typename I1, ..., typename IN>
     *      const_marray<T, ndim>(I1 len1, ..., IN lenN, unsigned align_at,
     *                            Layout layout = DEFAULT)
     *
     * 9.1) template <typename I1, ..., typename IN>
     *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const T& val,
     *                            unsigned align_at, Layout layout = DEFAULT)
     *
     * 9.1) template <typename I1, ..., typename IN>
     *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const unsigned_t& u,
     *                            unsigned align_at, Layout layout = DEFAULT)
     *
     *  The same as 6.1-6.3 and 7.1-7.3, except that the lengths are given as
     *  distinct parameters with possibly difference integral types. There must
     *  be exactly ndim length parameters.
     *
     * 10.1) template <typename U>
     *       const_marray<T, ndim>(const std::array<U, ndim>& len, const T* ptr,
     *                             Layout layout = DEFAULT)
     *
     * 10.2) template <typename U, typename V>
     *       const_marray<T, ndim>(const std::array<U, ndim>& len, const T* ptr,
     *                             const std::array<V, ndim>& stride)
     *
     *  Constructs an array with the given lengths and data layout, using ptr
     *  as the array data. The given pointer is referenced directly by the
     *  array, providing a "view" of an externally-allocated data segment.
     *  No internal allocation or copying is performed. In the second form the
     *  strides are specified directly, instead of inferred from the data
     *  layout. It is the user's reposnsibility to ensure that the information
     *  given results in legal accesses of ptr.
     *
     * 11.1) template <typename I1, ..., typename IN>
     *       const_marray<T, ndim>(I1 len1, ..., IN lenN, const T* ptr, Layout layout = DEFAULT)
     *
     * 11.2) template <typename I1, ..., typename IN, typename S1, ..., typename SN>
     *       const_marray<T, ndim>(I1 len1, ..., IN lenN, const T* ptr,
     *                             S1 stride1, ..., SN strideN)
     *
     *  The same as 10.1 and 10.2, except that the lengths and, in the second form
     *  the strides, are given as distinct parameters with possibly different
     *  integral types. There must be exactly ndim length and stride parameters.
     *
     * 12.1) template <...>
     *       const_marray<T, ndim>(array[][][]...)
     *
     * 12.2) template <...>
     *       const_marray<T, ndim>(array(...))
     *
     *  Creates a view from the given array reference as if it were repeatedly
     *  indexed with slice::all.
     *
     * Individual elements of the array may be accessed using one of two forms:
     *
     * 1) C-style access:
     *
     *  const_marray<double, 4> arr(4, 3, 9, 5);
     *  double el = arr[1][1][6][0];
     *
     * 2) FORTRAN-style access (using 0-based indices):
     *
     *  const_marray<double, 4> arr(4, 3, 9, 5);
     *  double el = arr(1, 1, 6, 0);
     *
     * These two methods are entirely equivalent. In addition, a (possibly
     * lower-dimensional) view of all or part of the array may be created
     * using this same notation. Instead of an integral index, a range_t
     * or slice::all_t argument denotes that the given range or all elements,
     * respectively, of a particular dimension are requested. For example:
     *
     *  using slice::all;
     *  const_marray<double, 4> arr(4, 3, 9, 5);
     *  const_marray<double, 2> sub_arr = arr[all][2][range(1,5)][1];
     *  // sub_arr has lengths (4, 4)
     *  double el = sub_arr[2][2];
     *
     * The last two constructors detailed above allow expressions of the type:
     *
     *  using slice::all;
     *  const_marray<double, 4> arr(4, 3, 9, 5);
     *  const_marray<double, 3> sub_arr = arr[all][2][range(1,5)];
     *  // this is identical to:
     *  // const_marray<double, 3> sub_arr = arr[all][2][range(1,5)][all];
     *  double el = sub_arr[2][2][0];
     *
     * Note that only the trailing dimensions may be implicitly sliced using
     * this notation.
     *
     * The type marray<T, ndim> provides non-const access to elements and for
     * the creation of non-const views.
     */
    template <typename T, unsigned ndim>
    class const_marray
    {
        static_assert(ndim > 0, "0-dimensional marrays are not allowed.");

        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

        public:
            typedef unsigned idx_type;
            typedef size_t size_type;
            typedef T value_type;
            typedef T* pointer;
            typedef const T* const_pointer;
            typedef T& reference;
            typedef const T& const_reference;

        protected:
            pointer data_ = NULL;
            size_type size_ = 0;
            std::array<idx_type,ndim> len_;
            std::array<size_type,ndim> stride_;
            bool is_view_ = false;
            Layout layout_ = DEFAULT;
            aligned_allocator<T, MARRAY_BASE_ALIGNMENT> alloc_;

            const_marray& operator=(const const_marray& other)
            {
                assert(other.len_ == len_);

                if (!data_) return *this;

                if (!other.is_view_ && !is_view_ && other.stride_ == stride_)
                {
                    std::copy(other.data_, other.data_+other.size_, data_);
                }
                else
                {
                    Iterator<idx_type, size_type> it(other.len_, other.stride_, stride_);
                    const_pointer a_ = other.data_;
                          pointer b_ =       data_;
                    while (it.nextIteration(a_, b_)) *b_ = *a_;
                }

                return *this;
            }

            /*
             * Default initialization.
             */

        public:
            const_marray() {}

        protected:
            void reset()
            {
                free();
                data_ = NULL;
                size_ = 0;
                len_.fill(0);
                stride_.fill(0);
                layout_ = DEFAULT;
                is_view_ = false;
            }

            MARRAY_TEST
            (
                marray<double, 3> a(2, 4, 5);
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                // check that repeated calls work
                a.reset();
                a.reset();
                assert(a.data() == NULL);
                assert(a.lengths() == make_array<idx_type>(0, 0, 0));
                assert(a.strides() == make_array<size_type>(0, 0, 0));
                assert(a.layout() == DEFAULT);
            )

            /*
             * Copy, move, and view.
             */

        public:
            const_marray(const const_marray& other)
            {
                reset(other);
            }

        protected:
            void reset(const const_marray& other)
            {
                free();

                size_ = other.size_;
                len_ = other.len_;
                stride_ = other.stride_;
                is_view_ = other.is_view_;
                layout_ = other.layout_;

                if (is_view_)
                {
                    data_ = other.data_;
                }
                else if (size_ > 0)
                {
                    data_ = alloc_.allocate(size_);
                    std::uninitialized_copy(other.data_, other.data_+size_, data_);
                }
            }

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 3> ca(a);
                const_marray<double, 3> va(a, construct_view);
                const_marray<double, 3> cva(va);

                assert(ca.data() != NULL);
                assert(ca.data() != a.data());
                assert(ca.lengths() == a.lengths());
                assert(ca.strides() == a.strides());
                assert(ca.layout() == a.layout());
                assert(!ca.isView());

                assert(cva.data() == a.data());
                assert(cva.lengths() == a.lengths());
                assert(cva.strides() == a.strides());
                assert(cva.layout() == a.layout());
                assert(cva.isView());
            )

        public:
            const_marray(const_marray&& other)
            {
                reset(std::move(other));
            }

        protected:
            void reset(const_marray&& other)
            {
                free();

                data_ = other.data_;
                size_ = other.size_;
                len_ = other.len_;
                stride_ = other.stride_;
                is_view_ = other.is_view_;
                layout_ = other.layout_;

                if (!is_view_)
                {
                    other.data_ = NULL;
                    other.reset();
                }
            }

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const double *data = a.data();
                const_marray<double, 3> ca(std::move(a));
                typedef const_marray<double, 3>::idx_type idx_type;
                typedef const_marray<double, 3>::size_type size_type;

                assert(ca.data() == data);
                assert(ca.lengths() == make_array<idx_type>(2, 4, 5));
                assert(ca.strides() == make_array<size_type>(1, 2, 2*4));
                assert(ca.layout() == COLUMN_MAJOR);
                assert(!ca.isView());

                assert(a.data() == NULL);
                assert(a.lengths() == make_array<idx_type>(0, 0, 0));
                assert(a.strides() == make_array<size_type>(0, 0, 0));
                assert(a.layout() == DEFAULT);
            )

        public:
            const_marray(const const_marray& other, const construct_view_t& cv)
            {
                reset(other, cv);
            }

        protected:
            void reset(const const_marray& other, const construct_view_t& cv)
            {
                free();

                data_ = other.data_;
                size_ = 0;
                len_ = other.len_;
                stride_ = other.stride_;
                is_view_ = true;
                layout_ = other.layout_;
            }

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 3> va(a, construct_view);

                assert(va.data() == a.data());
                assert(va.lengths() == a.lengths());
                assert(va.strides() == a.strides());
                assert(va.layout() == a.layout());
                assert(va.isView());
            )

        public:
            const_marray(const const_marray& other, const construct_copy_t& cp)
            {
                reset(other, cp);
            }

        protected:
            void reset(const const_marray& other, const construct_copy_t& cp)
            {
                if (other.is_view_)
                {
                    if (std::is_scalar<T>::value)
                    {
                        reset(other.len_, uninitialized, other.layout_);
                    }
                    else
                    {
                        reset(other.len_, T(), other.layout_);
                    }
                    *this = other;
                }
                else
                {
                    reset(other.len_, uninitialized, other.layout_);
                    std::uninitialized_copy(other.data_, other.data_+size_, data_);
                }
            }

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 3> va(a, construct_view);
                const_marray<double, 3> ca(va, construct_copy);

                assert(ca.data() != NULL);
                assert(ca.data() != a.data());
                assert(ca.lengths() == a.lengths());
                assert(ca.strides() == a.strides());
                assert(ca.layout() == a.layout());
                assert(!ca.isView());
            )

            /*
             * Direct initialization.
             */

        public:
            template <typename U>
            explicit const_marray(const std::array<U, ndim>& len, Layout layout = DEFAULT)
            {
                reset(len, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, ndim>& len, Layout layout = DEFAULT)
            {
                reset(len, T(), layout);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, ndim>& len, const T& val, Layout layout = DEFAULT)
            {
                reset(len, val, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, ndim>& len, const T& val, Layout layout = DEFAULT)
            {
                reset(len, uninitialized, layout);
                std::uninitialized_fill_n(data_, size_, val);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, ndim>& len, const uninitialized_t& u, Layout layout = DEFAULT)
            {
                reset(len, u, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, ndim>& len, const uninitialized_t& u, Layout layout = DEFAULT)
            {
                reset(len, u, 0, layout);
            }

        public:
            template <typename U>
            explicit const_marray(const std::array<U, ndim>& len, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, align_at, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, ndim>& len, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, T(), align_at, layout);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, ndim>& len, const T& val, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, val, align_at, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, ndim>& len, const T& val, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, uninitialized, align_at, layout);
                std::uninitialized_fill_n(data_, size_, val);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, ndim>& len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, u, align_at, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, ndim>& len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            {
                assert(align_at < ndim);

                free();

                std::copy_n(len.begin(), ndim, len_.begin());
                layout_ = layout;
                is_view_ = false;
                data_ = NULL;

                if (layout_ == ROW_MAJOR)
                {
                    stride_[ndim-1] = 1;
                    for (unsigned i = ndim-1;i > 0;i--)
                    {
                        if (i == align_at)
                        {
                            stride_[i-1] = detail::align(stride_[i]*len_[i], MARRAY_STRIDE_ALIGNMENT);
                        }
                        else
                        {
                            stride_[i-1] = stride_[i]*len_[i];
                        }
                    }
                    size_ = stride_[0]*len_[0];
                }
                else
                {
                    stride_[0] = 1;
                    for (unsigned i = 1;i < ndim;i++)
                    {
                        if (i == align_at)
                        {
                            stride_[i] = detail::align(stride_[i-1]*len_[i-1], MARRAY_STRIDE_ALIGNMENT);
                        }
                        else
                        {
                            stride_[i] = stride_[i-1]*len_[i-1];
                        }
                    }
                    size_ = stride_[ndim-1]*len_[ndim-1];
                }

                if (size_ > 0)
                {
                    data_ = alloc_.allocate(size_);
                }
            }

            MARRAY_TEST
            (
                marray<double, 3> a;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;

                a.reset(make_array(2, 4, 5));
                assert(a.data() != NULL);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(4*5, 5, 1));
                assert(a.layout() == DEFAULT);
                a.reset();
            )

            /*
             * Wrap external pointer.
             */

        public:
            template <typename U>
            const_marray(const std::array<U, ndim>& len, const_pointer ptr,
                         Layout layout = DEFAULT)
            {
                reset(len, const_cast<pointer>(ptr), layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, ndim>& len,
                       pointer ptr, Layout layout = DEFAULT)
            {
                if (layout == ROW_MAJOR)
                {
                    stride_[ndim-1] = 1;
                    for (unsigned i = ndim-1;i > 0;i--)
                    {
                        stride_[i-1] = stride_[i]*len[i];
                    }
                }
                else
                {
                    stride_[0] = 1;
                    for (unsigned i = 1;i < ndim;i++)
                    {
                        stride_[i] = stride_[i-1]*len[i-1];
                    }
                }

                reset(len, ptr, stride_);
                layout_ = layout;
            }

        public:
            template <typename U, typename V>
            const_marray(const std::array<U, ndim>& len, const_pointer ptr,
                         const std::array<V, ndim>& stride)
            {
                reset(len, const_cast<pointer>(ptr), stride);
            }

        protected:
            template <typename U, typename V>
            void reset(const std::array<U, ndim>& len, pointer ptr,
                       const std::array<V, ndim>& stride)
            {
                free();

                std::copy_n(len.begin(), ndim, len_.begin());
                std::copy_n(stride.begin(), ndim, stride_.begin());
                data_ = ptr;
                size_ = 0;
                layout_ = DEFAULT;
                is_view_ = true;
            }

            MARRAY_TEST
            (
                double p;
                marray<double, 3> a;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;

                a.reset(make_array(2, 4, 5), &p);
                assert(a.data() == &p);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(4*5, 5, 1));
                assert(a.layout() == DEFAULT);
                a.reset();

                a.reset(make_array(2, 4, 5), &p, make_array(25, 6, 1));
                assert(a.data() == &p);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(25, 6, 1));
                assert(a.layout() == DEFAULT);
                a.reset();
            )

            /*
             * Polymorphic initialization.
             */

        public:
            template <typename... Args>
            explicit const_marray(typename std::enable_if<detail::are_reset_args<T, ndim, idx_type, Args...>::value,idx_type>::type len0, Args&&... args)
            {
                detail::are_reset_args<T, ndim, idx_type, Args...>(*this, len0, std::forward<Args>(args)...);
            }

        protected:
            template <typename... Args>
            void reset(typename std::enable_if<detail::are_reset_args<T, ndim, idx_type, Args...>::value,idx_type>::type len0, Args&&... args)
            {
                detail::are_reset_args<T, ndim, idx_type, Args...>(*this, len0, std::forward<Args>(args)...);
            }

            MARRAY_TEST
            (
                double p;
                marray<double, 3> a;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;

                a.reset(2, 4, 5);
                assert(a.data() != NULL);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(4*5, 5, 1));
                assert(a.layout() == DEFAULT);
                a.reset();

                a.reset(2, 4, 5, &p);
                assert(a.data() == &p);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(4*5, 5, 1));
                assert(a.layout() == DEFAULT);
                a.reset();

                a.reset(2, 4, 5, &p, 25, 6, 1);
                assert(a.data() == &p);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(25, 6, 1));
                assert(a.layout() == DEFAULT);
                a.reset();
            )

            /*
             * Construct view from ref or slice.
             */

        public:
            template <unsigned ndim_>
            const_marray(const const_marray_ref<T, ndim_, ndim_-ndim+1>& other)
            : data_(other.array_.data_+other.idx), is_view_(true), layout_(other.array_.layout_)
            {
                std::copy_n(other.array_.len_.begin()+ndim_-ndim, ndim, len_.begin());
                std::copy_n(other.array_.stride_.begin()+ndim_-ndim, ndim, stride_.begin());
            }

            MARRAY_TEST
            (
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 2> b = a[0];

                assert(a.data() == b.data());
                assert(b.lengths() == make_array<idx_type>(4, 5));
                assert(b.strides() == make_array<size_type>(2, 2*4));
                assert(b.layout() == COLUMN_MAJOR);
            )

            template <unsigned ndim_, unsigned newdim_>
            const_marray(const const_marray_slice<T, ndim_, ndim_-ndim+1+newdim_, newdim_>& other)
            : data_(other.array_.data_+other.idx), is_view_(true), layout_(other.array_.layout_)
            {
                for (unsigned i = 0;i < newdim_;i++)
                {
                    len_[i] = other.lens[i];
                    stride_[i] = other.array_.stride_[other.dims[i]];
                }
                std::copy_n(other.array_.len_.begin()+ndim_-ndim+newdim_, ndim-newdim_, len_.begin()+newdim_);
                std::copy_n(other.array_.stride_.begin()+ndim_-ndim+newdim_, ndim-newdim_, stride_.begin()+newdim_);
            }

            MARRAY_TEST
            (
                using slice::all;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                const_marray<double, 4> a(2, 1, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 3> b = a[0][all];

                assert(a.data() == b.data());
                assert(b.lengths() == make_array<idx_type>(1, 4, 5));
                assert(b.strides() == make_array<size_type>(2, 2*1, 2*1*4));
                assert(b.layout() == COLUMN_MAJOR);
            )

            /*
             * Destruction.
             */

        public:
            ~const_marray()
            {
                free();
            }

        protected:
            void free()
            {
                if (!is_view_ && data_)
                {
                    for (size_type i = 0;i < size_;i++)
                    {
                        std::allocator_traits<decltype(alloc_)>::destroy(alloc_, data_+i);
                    }
                    alloc_.deallocate(data_, size_);
                }
            }

            /*
             * Views.
             */

        public:
            bool isView() const
            {
                return is_view_;
            }

            const_marray<T, ndim> view() const
            {
                return const_marray<T, ndim>(*this, construct_view);
            }

            template <typename U>
            const_marray<T, ndim> permute(const std::array<U, ndim>& perm) const
            {
                const_marray<T, ndim> view;
                view.data_ = data_;
                view.is_view_ = true;
                view.layout_ = layout_;

                for (unsigned i = 0;i < ndim;i++)
                {
                    assert(0 <= perm[i] && perm[i] < ndim);
                    for (unsigned j = 0;j < i;j++) assert(perm[i] != perm[j]);
                }

                for (unsigned i = 0;i < ndim;i++)
                {
                    view.len_[i] = len_[perm[i]];
                    view.stride_[i] = stride_[perm[i]];
                }

                return view;
            }

            template <typename... Args>
            typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, const_marray<T, ndim>>::type
            permute(Args&&... args) const
            {
                return permute(make_array(std::forward<Args>(args)...));
            }

            /*
             * Access view.
             */

            const_marray<T, ndim-1> front(unsigned dim) const
            {
                assert(dim < ndim);
                assert(len_[dim] > 0);

                const_marray<T, ndim-1> view;
                view.data_ = data_;
                copy_n(stride_.begin(), dim, view.stride_.begin());
                copy_n(stride_.begin()+dim+1, ndim-dim-1, view.stride_.begin()+dim);
                copy_n(len_.begin(), dim, view.len_.begin());
                copy_n(len_.begin()+dim+1, ndim-dim-1, view.len_.begin()+dim);
                view.is_view_ = true;
                view.layout_ = layout_;

                return view;
            }

            const_marray<T, ndim-1> back(unsigned dim) const
            {
                const_marray<T, ndim-1> view = front(dim);
                view.data_ += (len_[dim]-1)*stride_[dim];
                return view;
            }

            /*
             * Access sub-arrays.
             */

            const_marray_ref<T, ndim, 2> operator[](idx_type i) const
            {
                assert(i < len_[0]);
                return const_marray_ref<T, ndim, 2>(*this, (size_type)0, i);
            }

            template <typename I>
            const_marray_slice<T, ndim, 2, 1> operator[](const range_t<I>& x) const
            {
                assert(x.front() >= 0 && x.back() <= len_[0]);
                return const_marray_slice<T, ndim, 2, 1>(*this, (size_type)0, {}, {}, x);
            }

            const_marray_slice<T, ndim, 2, 1> operator[](const slice::all_t& x) const
            {
                return const_marray_slice<T, ndim, 2, 1>(*this, (size_type)0, {}, {}, range(idx_type(), len_[0]));
            }

            /*
             * Access elements and sub-arrays.
             */

            template <typename... Indices>
            typename std::enable_if<sizeof...(Indices) == ndim, typename detail::return_type<true, T, Indices...>::type>::type
            operator()(const Indices&... idx) const
            {
                return detail::get_slice<typename detail::return_type<true, T, Indices...>::type, Indices...>()(*this, idx...);
            }

            /*
             * Member accessors.
             */

            const_pointer data() const
            {
                return data_;
            }

            const idx_type& length(unsigned dim) const
            {
                return len_[dim];
            }

            const std::array<idx_type, ndim>& lengths() const
            {
                return len_;
            }

            const size_type& stride(unsigned dim) const
            {
                return stride_[dim];
            }

            const std::array<size_type, ndim>& strides() const
            {
                return stride_;
            }

            Layout layout() const
            {
                return layout_;
            }
    };

    /*
     * This class is identical to const_marray except that modification of the
     * array elements is allowed. Objects of this type may be implicitly
     * converted to const_marray, but not vice-versa. This means that, for example
     * a function which takes a const array reference must have a parameter of the
     * type 'const const_marray<T, ndim>&' instead of 'const marray<T, ndim>&'.
     */
    template <typename T, unsigned ndim>
    class marray : public const_marray<T, ndim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, typename U, typename... Args> friend struct detail::marray_reset;
        template <typename T_, unsigned ndim_, unsigned dim_, typename U, typename... Args> friend struct detail::marray_resize;

        public:
            using typename const_marray<T, ndim>::idx_type;
            using typename const_marray<T, ndim>::size_type;
            using typename const_marray<T, ndim>::value_type;
            using typename const_marray<T, ndim>::pointer;
            using typename const_marray<T, ndim>::const_pointer;
            using typename const_marray<T, ndim>::reference;
            using typename const_marray<T, ndim>::const_reference;

        protected:
            using const_marray<T, ndim>::data_;
            using const_marray<T, ndim>::size_;
            using const_marray<T, ndim>::len_;
            using const_marray<T, ndim>::stride_;
            using const_marray<T, ndim>::is_view_;
            using const_marray<T, ndim>::layout_;

            /*
             * Cheater functions to allow non-const versions of view-creating
             * methods to use return values from const versions.
             */
            marray(const_marray<T, ndim>&& other)
            : const_marray<T, ndim>(std::move(other)) {}

            marray(const const_marray<T, ndim>& other)
            : const_marray<T, ndim>(other) {}

        public:
            /*
             * Default initialization.
             */

            marray() {}

            /*
             * Copy, move, and view.
             */

            marray(marray& other)
            : const_marray<T, ndim>(other) {}

            marray(marray&& other)
            : const_marray<T, ndim>(std::move(other)) {}

            marray(marray& other, const construct_view_t& cv)
            : const_marray<T, ndim>(other, cv) {}

            marray(const marray& other, const construct_copy_t& cp)
            : const_marray<T, ndim>(other, cp) {}

            /*
             * Direct initialization.
             */

            template <typename U>
            explicit marray(const std::array<U, ndim>& len, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, layout) {}

            template <typename U>
            marray(const std::array<U, ndim>& len, const T& val, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, val, layout) {}

            template <typename U>
            marray(const std::array<U, ndim>& len, const uninitialized_t& u, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, u, layout) {}

            template <typename U>
            marray(const std::array<U, ndim>& len, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, align_at, layout) {}

            template <typename U>
            marray(const std::array<U, ndim>& len, const T& val, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, val, align_at, layout) {}

            template <typename U>
            marray(const std::array<U, ndim>& len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, u, align_at, layout) {}

            /*
             * Wrap external pointer.
             */

            template <typename U>
            marray(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, ptr, layout) {}

            template <typename U, typename V>
            marray(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
            : const_marray<T, ndim>(len, ptr, stride) {}

            /*
             * Polymorphic initialization.
             */

            template <typename... Args>
            marray(typename std::enable_if<detail::are_reset_args<T, ndim, idx_type, Args...>::value,idx_type>::type len0, Args&&... args)
            : const_marray<T, ndim>(len0, std::forward<Args>(args)...) {}

            /*
             * Construct view from ref or slice.
             */

            template <unsigned ndim_>
            marray(const marray_ref<T, ndim_, ndim_-ndim+1>& other)
            : const_marray<T, ndim>(other) {}

            template <unsigned ndim_, unsigned newdim_>
            marray(const marray_slice<T, ndim_, ndim_-ndim+1+newdim_, newdim_>& other)
            : const_marray<T, ndim>(other) {}

            /*
             * Reset and clear.
             */

            using const_marray<T, ndim>::reset;

            void clear()
            {
                reset();
            }

            /*
             * Resize.
             */

            template <typename U>
            void resize(const std::array<U, ndim>& len, const T& val = T())
            {
                const_marray<T, ndim> a(std::move(*this));
                reset(len, val, layout_);
                const_marray<T, ndim> b = view();

                if (a.data_ && b.data_)
                {
                    /*
                     * It is OK to change the geometry of 'a' even if it is not
                     * a view since it is about to go out of scope.
                     */
                    for (unsigned i = 0;i < ndim;i++)
                    {
                        a.len_[i] = b.len_[i] = std::min(a.len_[i], b.len_[i]);
                    }
                    b = a;
                }
            }

            template <typename... Args>
            void resize(typename std::enable_if<detail::are_resize_args<T, ndim, idx_type, Args...>::value,idx_type>::type len0, Args&&... args)
            {
                detail::are_resize_args<T, ndim, idx_type, Args...>(*this, len0, std::forward<Args>(args)...);
            }

            /*
             * View -> copy.
             */

            void unView()
            {
                if (!is_view_) return;

                const_marray<T, ndim> old(*this);
                reset(len_, layout_);
                *this = old;
            }

            /*
             * Views.
             */

            using const_marray<T, ndim>::view;

            marray<T, ndim> view()
            {
                return marray<T, ndim>(*this, construct_view);
            }

            using const_marray<T, ndim>::permute;

            template <typename U>
            marray<T, ndim> permute(const std::array<U, ndim>& perm)
            {
                return const_marray<T, ndim>::permute(perm);
            }

            template <typename... Args>
            typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, marray<T, ndim>>::type
            permute(Args&&... args)
            {
                return const_marray<T, ndim>::permute(std::forward<Args>(args)...);
            }

            /*
             * Modify.
             */

            void push_back(unsigned dim, const const_marray<T, ndim-1>& x)
            {
                assert(dim < ndim);

                for (unsigned i = 0, j = 0;i < ndim;i++)
                {
                    if (i != dim)
                    {
                        assert(len_[i] == x.len_[j++]);
                    }
                }

                std::array<idx_type, ndim> len = len_;
                len[dim]++;
                resize(len);
                back() = x;
            }

            void pop_back(unsigned dim)
            {
                assert(dim < ndim);
                assert(len_[dim] > 0);

                std::array<idx_type, ndim> len = len_;
                len[dim]--;
                resize(len);
            }

            /*
             * Access view.
             */

            using const_marray<T, ndim>::front;

            marray<T, ndim-1> front(unsigned dim)
            {
                return const_marray<T, ndim>::front(dim);
            }

            using const_marray<T, ndim>::back;

            marray<T, ndim-1> back(unsigned dim)
            {
                return const_marray<T, ndim>::back(dim);
            }

            /*
             * Assign.
             */

            marray& operator=(const marray& other)
            {
                const_marray<T, ndim>::operator=(other);
                return *this;
            }

            marray& operator=(const const_marray<T, ndim>& other)
            {
                const_marray<T, ndim>::operator=(other);
                return *this;
            }

            marray& operator=(const T& x)
            {
                if (is_view_)
                {
                    Iterator<idx_type, size_type> it(len_, stride_);
                    pointer p = data_;
                    while (it.nextIteration(p)) *p = x;
                }
                else
                {
                    std::fill(data_, data_+size_, x);
                }
                return *this;
            }

            /*
             * Access sub-arrays.
             */

            using const_marray<T, ndim>::operator[];

            marray_ref<T, ndim, 2> operator[](idx_type i)
            {
                assert(i < len_[0]);
                return marray_ref<T, ndim, 2>(*this, (size_type)0, i);
            }

            template <typename I>
            marray_slice<T, ndim, 2, 1> operator[](const range_t<I>& x)
            {
                assert(x.front() >= 0 && x.back() <= len_[0]);
                return marray_slice<T, ndim, 2, 1>(*this, (size_type)0, {}, {}, x);
            }

            marray_slice<T, ndim, 2, 1> operator[](const slice::all_t& x)
            {
                return marray_slice<T, ndim, 2, 1>(*this, (size_type)0, {}, {}, range(idx_type(), len_[0]));
            }

            /*
             * Access elements and sub-arrays.
             */

            using const_marray<T, ndim>::operator();

            template <typename... Indices>
            typename std::enable_if<sizeof...(Indices) == ndim, typename detail::return_type<false, T, Indices...>::type>::type
            operator()(const Indices&... idx)
            {
                return detail::get_slice<typename detail::return_type<false, T, Indices...>::type, Indices...>()(*this, idx...);
            }

            /*
             * Member accessors.
             */

            using const_marray<T, ndim>::data;

            pointer data()
            {
                return data_;
            }

            /*
             * Swap.
             */

            void swap(marray&& other)
            {
                swap(other);
            }

            void swap(marray& other)
            {
                using std::swap;
                swap(data_,      other.data_);
                swap(size_,      other.size_);
                swap(len_,       other.len_);
                swap(stride_,    other.stride_);
                swap(is_view_,   other.is_view_);
                swap(layout_,    other.layout_);
            }

            friend void swap(marray&& a, marray&& b)
            {
                a.swap(b);
            }

            friend void swap(marray& a, marray& b)
            {
                a.swap(b);
            }
    };

    /*
     * This class represents a one-dimensional array of elements of type T.
     * This specialization is provided to optimize the implementation of
     * linear arrays.
     */
    template <typename T>
    class const_marray<T, 1>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            typedef size_t idx_type;
            typedef size_t size_type;
            typedef T value_type;
            typedef T* pointer;
            typedef const T* const_pointer;
            typedef T& reference;
            typedef const T& const_reference;

        protected:
            pointer data_ = NULL;
            idx_type size_ = 0;
            std::array<idx_type, 1> len_ = {{0}};
            std::array<size_type, 1> stride_ = {{1}};
            bool is_view_ = false;
            Layout layout_ = DEFAULT;
            aligned_allocator<T, MARRAY_BASE_ALIGNMENT> alloc_;

            const_marray& operator=(const const_marray& other)
            {
                assert(length() == other.length());

                if (!other.is_view_ && !is_view_)
                {
                    std::copy_n(other.data_, length(), data_);
                }
                else
                {
                    const_pointer a_ = other.data_;
                          pointer b_ =       data_;

                    size_type as =       stride();
                    size_type bs = other.stride();

                    for (size_type i = 0;i < length();i++)
                    {
                        b_[i*bs] = a_[i*as];
                    }
                }

                return *this;
            }

            /*
             * Default initialization.
             */

        public:
            const_marray() {}

        protected:
            void reset()
            {
                free();
                data_ = NULL;
                size_ = 0;
                len_.fill(0);
                stride_.fill(0);
                is_view_ = false;
            }

            /*
             * Copy, move, and view.
             */

        public:
            const_marray(const const_marray& other)
            {
                reset(other);
            }

        protected:
            void reset(const const_marray& other)
            {
                free();

                size_ = other.size_;
                len_ = other.len_;
                stride_ = other.stride_;
                is_view_ = other.is_view_;

                if (is_view_)
                {
                    data_ = other.data_;
                }
                else if (size_ > 0)
                {
                    data_ = alloc_.allocate(size_);
                    std::uninitialized_copy(other.data_, other.data_+size_, data_);
                }
            }

        public:
            const_marray(const_marray&& other)
            {
                reset(std::move(other));
            }

        protected:
            void reset(const_marray&& other)
            {
                free();

                data_ = other.data_;
                size_ = other.size_;
                len_ = other.len_;
                stride_ = other.stride_;
                is_view_ = other.is_view_;

                if (!is_view_)
                {
                    other.data_ = NULL;
                    other.reset();
                }
            }

        public:
            const_marray(const const_marray& other, const construct_view_t& cv)
            {
                reset(other, cv);
            }

        protected:
            void reset(const const_marray& other, const construct_view_t& cv)
            {
                free();

                data_ = other.data_;
                size_ = 0;
                len_ = other.len_;
                stride_ = other.stride_;
                is_view_ = true;
            }

        public:
            const_marray(const const_marray& other, const construct_copy_t& cp)
            {
                reset(other, cp);
            }

        protected:
            void reset(const const_marray& other, const construct_copy_t& cp)
            {
                if (other.is_view_)
                {
                    if (std::is_scalar<T>::value)
                    {
                        reset(other.len_, uninitialized);
                    }
                    else
                    {
                        reset(other.len_, T());
                    }
                    *this = other;
                }
                else
                {
                    reset(other.len_, uninitialized);
                    std::uninitialized_copy(other.data_, other.data_+size_, data_);
                }
            }

            /*
             * Direct initialization.
             */

        public:
            template <typename U>
            explicit const_marray(const std::array<U, 1>& len, Layout layout = DEFAULT)
            {
                reset(len, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, 1>& len, Layout layout = DEFAULT)
            {
                reset(len, T(), layout);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, 1>& len, const T& val, Layout layout = DEFAULT)
            {
                reset(len, val, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, 1>& len, const T& val, Layout layout = DEFAULT)
            {
                reset(len, uninitialized, layout);
                std::uninitialized_fill_n(data_, size_, val);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, 1>& len, const uninitialized_t& u, Layout layout = DEFAULT)
            {
                reset(len, u, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, 1>& len, const uninitialized_t& u, Layout layout = DEFAULT)
            {
                reset(len, u, 0, layout);
            }

        public:
            template <typename U>
            explicit const_marray(const std::array<U, 1>& len, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, align_at, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, 1>& len, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, T(), align_at, layout);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, 1>& len, const T& val, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, val, align_at, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, 1>& len, const T& val, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, uninitialized, align_at, layout);
                std::uninitialized_fill_n(data_, size_, val);
            }

        public:
            template <typename U>
            const_marray(const std::array<U, 1>& len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len, u, align_at, layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, 1>& len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(len[0], u, align_at, layout);
            }

        public:
            explicit const_marray(idx_type n, Layout layout = DEFAULT)
            {
                reset(n, layout);
            }

        protected:
            void reset(idx_type n, Layout layout = DEFAULT)
            {
                reset(n, T(), layout);
            }

        public:
            const_marray(idx_type n, const T& val, Layout layout = DEFAULT)
            {
                reset(n, val, layout);
            }

        protected:
            void reset(idx_type n, const T& val, Layout layout = DEFAULT)
            {
                reset(n, uninitialized, layout);
                std::uninitialized_fill_n(data_, size_, val);
            }

        public:
            const_marray(idx_type n, const uninitialized_t& u, Layout layout = DEFAULT)
            {
                reset(n, u, layout);
            }

        protected:
            void reset(idx_type n, const uninitialized_t& u, Layout layout = DEFAULT)
            {
                reset(n, u, 0, layout);
            }

        public:
            explicit const_marray(idx_type n, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(n, align_at, layout);
            }

        protected:
            void reset(idx_type n, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(n, T(), align_at, layout);
            }

        public:
            const_marray(idx_type n, const T& val, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(n, val, align_at, layout);
            }

        protected:
            void reset(idx_type n, const T& val, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(n, uninitialized, align_at, layout);
                std::uninitialized_fill_n(data_, size_, val);
            }

        public:
            const_marray(idx_type n, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            {
                reset(n, u, align_at, layout);
            }

        protected:
            void reset(idx_type n, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            {
                assert(align_at == 0);

                free();

                size_ = n;
                len_[0] = n;
                stride_[0] = 1;
                is_view_ = false;
                data_ = NULL;

                if (size_ > 0)
                {
                    data_ = alloc_.allocate(size_);
                }
            }

            /*
             * Wrap external pointer.
             */

        public:
            template <typename U>
            const_marray(const std::array<U, 1>& len, const_pointer ptr,
                         Layout layout = DEFAULT)
            {
                reset(len, const_cast<pointer>(ptr), layout);
            }

        protected:
            template <typename U>
            void reset(const std::array<U, 1>& len,
                       pointer ptr, Layout layout = DEFAULT)
            {
                reset(len[0], ptr, 1);
            }

            const_marray& operator=(const const_marray& other)
            {
                assert(size_ == other.size_);
                const_pointer a_ = other.data_;
                      pointer b_ =       data_;
                for (size_type i = 0;i < size_;i++)
                {
                    *b_ = *a_;
                    a_ += other.stride_;
                    b_ +=       stride_;
                }
                return *this;
            }

        public:
            template <typename U, typename V>
            const_marray(const std::array<U, 1>& len, const_pointer ptr,
                         const std::array<V, 1>& stride)
            {
                reset(len, const_cast<pointer>(ptr), stride);
            }

        protected:
            template <typename U, typename V>
            void reset(const std::array<U, 1>& len, pointer ptr,
                       const std::array<V, 1>& stride)
            {
                reset(len[0], ptr, stride[0]);
            }

        public:
            template <typename U>
            const_marray(idx_type len, const_pointer ptr, Layout layout = DEFAULT)
            {
                reset(len, const_cast<pointer>(ptr), layout);
            }

        protected:
            template <typename U>
            void reset(idx_type len, pointer ptr, Layout layout = DEFAULT)
            {
                reset(len, ptr, 1);
            }

        public:
            template <typename U, typename V>
            const_marray(idx_type len, pointer ptr, size_type stride)
            {
                reset(len, const_cast<pointer>(ptr), stride);
            }

        protected:
            template <typename U, typename V>
            void reset(idx_type len, pointer ptr, size_type stride)
            {
                free();

                len_ = len;
                stride_ = stride;
                data_ = ptr;
                size_ = 0;
                is_view_ = true;
            }

            /*
             * Construct view from ref or slice.
             */

        public:
            template <unsigned ndim_>
            const_marray(const const_marray_ref<T, ndim_, ndim_>& other)
            : data_(other.array_.data_+other.idx), len_({other.array_.len_[ndim_-1]}),
              stride_({other.array_.stride_[ndim_-1]}), is_view_(true) {}

            /*
             * Destruction.
             */

        public:
            ~const_marray()
            {
                free();
            }

        protected:
            void free()
            {
                if (!is_view_ && data_)
                {
                    for (size_type i = 0;i < size_;i++)
                    {
                        std::allocator_traits<decltype(alloc_)>::destroy(alloc_, data_+i);
                    }
                    alloc_.deallocate(data_, size_);
                }
            }

            /*
             * Views.
             */

        public:
            bool isView() const
            {
                return is_view_;
            }

            const_marray<T, 1> view() const
            {
                return const_marray<T, 1>(*this, construct_view);
            }

            template <typename U>
            const_marray<T, 1> permute(const std::array<U, 1>& perm) const
            {
                assert(perm[0] == 0);
                return view();
            }

            const_marray<T, 1> permute(unsigned perm) const
            {
                assert(perm == 0);
                return view();
            }

            /*
             * Access view.
             */

            const_reference front(unsigned dim) const
            {
                assert(dim == 0);
                return front();
            }

            const_reference back(unsigned dim) const
            {
                assert(dim == 0);
                return back();
            }

            const_reference front() const
            {
                assert(length() > 0);
                return data_[0];
            }

            const_reference back() const
            {
                assert(length() > 0);
                return data_[(length()-1)*stride()];
            }

            /*
             * Access elements and sub-arrays.
             */

            const_reference operator[](idx_type i) const
            {
                assert(i < length());
                return data_[i*stride()];
            }

            template <typename I>
            const_marray<T, 1> operator[](const range_t<I>& x) const
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= length());
                return const_marray<T, 1>(x.size(), data_+x.front()*stride(), stride());
            }

            const const_marray<T, 1>& operator[](const slice::all_t& x) const
            {
                return *this;
            }

            const_reference operator()(idx_type i) const
            {
                return (*this)[i];
            }

            template <typename I>
            const_marray<T, 1> operator()(const range_t<I>& x) const
            {
                return (*this)[x];
            }

            const const_marray<T, 1>& operator()(const slice::all_t& x) const
            {
                return *this;
            }

            /*
             * Member accessors.
             */

            const_pointer data() const
            {
                return data_;
            }

            idx_type length() const
            {
                return len_[0];
            }

            idx_type length(unsigned dim) const
            {
                return len_[0];
            }

            std::array<idx_type, 1> lengths() const
            {
                return len_;
            }

            size_type stride() const
            {
                return stride_[0];
            }

            size_type stride(unsigned dim) const
            {
                return stride_[0];
            }

            std::array<size_type, 1> strides() const
            {
                return stride_;
            }

            /*
             * Provided for similarity to std::vector.
             */
            size_type size() const
            {
                return len_[0];
            }
    };

    /*
     * This class is identical to const_marray<T, 1> except that it allows
     * modification of the array elements.
     */
    template <typename T>
    class marray<T, 1> : public const_marray<T, 1>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            using typename const_marray<T, 1>::idx_type;
            using typename const_marray<T, 1>::size_type;
            using typename const_marray<T, 1>::value_type;
            using typename const_marray<T, 1>::pointer;
            using typename const_marray<T, 1>::const_pointer;
            using typename const_marray<T, 1>::reference;
            using typename const_marray<T, 1>::const_reference;

        protected:
            using const_marray<T, 1>::data_;
            using const_marray<T, 1>::size_;
            using const_marray<T, 1>::len_;
            using const_marray<T, 1>::stride_;
            using const_marray<T, 1>::is_view_;

            /*
             * Cheater functions to allow non-const versions of view-creating
             * methods to use return values from const versions.
             */
            marray(const_marray<T, 1>&& other)
            : const_marray<T, 1>(std::move(other)) {}

            marray(const const_marray<T, 1>& other)
            : const_marray<T, 1>(other) {}

        public:
            using const_marray<T, 1>::length;
            using const_marray<T, 1>::stride;

            /*
             * Default initialization.
             */

            marray() {}

            /*
             * Copy, move, and view.
             */

            marray(marray& other)
            : const_marray<T, 1>(other) {}

            marray(marray&& other)
            : const_marray<T, 1>(std::move(other)) {}

            marray(marray& other, const construct_view_t& cv)
            : const_marray<T, 1>(other, cv) {}

            marray(const marray& other, const construct_copy_t& cp)
            : const_marray<T, 1>(other, cp) {}

            /*
             * Direct initialization.
             */

            template <typename U>
            explicit marray(const std::array<U, 1>& len, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, layout) {}

            template <typename U>
            marray(const std::array<U, 1>& len, const T& val, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, val, layout) {}

            template <typename U>
            marray(const std::array<U, 1>& len, const uninitialized_t& u, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, u, layout) {}

            template <typename U>
            marray(const std::array<U, 1>& len, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, align_at, layout) {}

            template <typename U>
            marray(const std::array<U, 1>& len, const T& val, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, val, align_at, layout) {}

            template <typename U>
            marray(const std::array<U, 1>& len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, u, align_at, layout) {}

            explicit marray(idx_type len, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, layout) {}

            marray(idx_type len, const T& val, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, val, layout) {}

            marray(idx_type len, const uninitialized_t& u, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, u, layout) {}

            marray(idx_type len, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, align_at, layout) {}

            marray(idx_type len, const T& val, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, val, align_at, layout) {}

            marray(idx_type len, const uninitialized_t& u, unsigned align_at, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, u, align_at, layout) {}

            /*
             * Wrap external pointer.
             */

            template <typename U>
            marray(const std::array<U, 1>& len, pointer ptr, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, ptr, layout) {}

            template <typename U, typename V>
            marray(const std::array<U, 1>& len, pointer ptr, const std::array<V, 1>& stride)
            : const_marray<T, 1>(len, ptr, stride) {}

            template <typename U>
            marray(idx_type len, pointer ptr, Layout layout = DEFAULT)
            : const_marray<T, 1>(len, ptr, layout) {}

            template <typename U, typename V>
            marray(idx_type len, pointer ptr, size_type stride)
            : const_marray<T, 1>(len, ptr, stride) {}

            /*
             * Construct view from ref or slice.
             */

            template <unsigned ndim_>
            marray(const marray_ref<T, ndim_, ndim_>& other)
            : const_marray<T, 1>(other) {}

            /*
             * Reset and clear.
             */

            using const_marray<T, 1>::reset;

            void clear()
            {
                reset();
            }

            /*
             * Resize.
             */

            template <typename U>
            void resize(const std::array<U, 1>& len, const T& val = T())
            {
                resize(len[0], val);
            }

            void resize(idx_type n, const T& val = T())
            {
                if (n <= length())
                {
                    len_[0] = n;
                    if (!is_view_) size_ = n;
                    return;
                }

                const_marray<T, 1> a(std::move(*this));
                reset(n, val);
                const_marray<T, 1> b = view();

                if (a.data_ && b.data_)
                {
                    /*
                     * It is OK to change the geometry of 'a' even if it is not
                     * a view since it is about to go out of scope.
                     */
                    a.len_[0] = b.len_[0] = std::min(a.len_[0], b.len_[0]);
                    b = a;
                }
                assert(data_ || length() == 0);
            }

            /*
             * View -> copy.
             */

            void unView()
            {
                if (!is_view_) return;

                const_marray<T, 1> old(*this);
                reset(len_);
                *this = old;
            }

            /*
             * Views.
             */

            using const_marray<T, 1>::view;

            marray<T, 1> view()
            {
                return marray<T, 1>(*this, construct_view);
            }

            using const_marray<T, 1>::permute;

            template <typename U>
            marray<T, 1> permute(const std::array<U, 1>& perm)
            {
                return const_marray<T, 1>::permute(perm);
            }

            marray<T, 1> permute(unsigned perm)
            {
                return const_marray<T, 1>::permute(perm);
            }

            /*
             * Modify.
             */

            void push_back(unsigned dim, const T& x)
            {
                assert(dim == 0);
                push_back(x);
            }

            void push_back(const T& x)
            {
                resize(length()+1);
                back() = x;
            }

            void pop_back(unsigned dim)
            {
                assert(dim == 0);
                pop_back();
            }

            void pop_back()
            {
                resize(length()-1);
            }

            /*
             * Access view.
             */

            using const_marray<T, 1>::front;

            reference front(unsigned dim)
            {
                return const_cast<reference>(const_marray<T, 1>::front(dim));
            }

            reference front()
            {
                return const_cast<reference>(const_marray<T, 1>::front());
            }

            using const_marray<T, 1>::back;

            reference back(unsigned dim)
            {
                return const_cast<reference>(const_marray<T, 1>::back(dim));
            }

            reference back()
            {
                return const_cast<reference>(const_marray<T, 1>::back());
            }

            /*
             * Assign.
             */

            marray& operator=(const marray& other)
            {
                const_marray<T, 1>::operator=(other);
                return *this;
            }

            marray& operator=(const const_marray<T, 1>& other)
            {
                const_marray<T, 1>::operator=(other);
                return *this;
            }

            marray& operator=(const T& x)
            {
                if (!is_view_)
                {
                    std::fill_n(data_, length(), x);
                }
                else
                {
                    for (size_type i = 0;i < length();i++)
                    {
                        data_[i*stride()] = x;
                    }
                }
                return *this;
            }

            /*
             * Access elements and sub-arrays.
             */

            using const_marray<T, 1>::operator[];

            reference operator[](idx_type i)
            {
                assert(i < length());
                return data_[i*stride()];
            }

            template <typename I>
            marray<T, 1> operator[](const range_t<I>& x)
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= size_);
                return marray<T, 1>(x.size(), data_+x.front()*stride_[0], stride_[0]);
            }

            marray<T, 1>& operator[](const slice::all_t& x)
            {
                return *this;
            }

            using const_marray<T, 1>::operator();

            reference operator()(idx_type i)
            {
                return (*this)[i];
            }

            template <typename I>
            marray<T, 1> operator()(const range_t<I>& x)
            {
                return (*this)[x];
            }

            marray<T, 1>& operator()(const slice::all_t& x)
            {
                return *this;
            }

            /*
             * Member accessors.
             */

            using const_marray<T, 1>::data;

            pointer data()
            {
                return data_;
            }

            /*
             * Swap.
             */

            void swap(marray&& other)
            {
                swap(other);
            }

            void swap(marray& other)
            {
                using std::swap;
                swap(data_,      other.data_);
                swap(size_,      other.size_);
                swap(len_,       other.len_);
                swap(stride_,    other.stride_);
                swap(is_view_,   other.is_view_);
            }

            friend void swap(marray&& a, marray&& b)
            {
                a.swap(b);
            }

            friend void swap(marray& a, marray& b)
            {
                a.swap(b);
            }
    };

    /*
     * Convenient names for 1- and 2-dimensional array types.
     */
    template <typename T> using const_row = const_marray<T, 1>;
    template <typename T> using row = marray<T, 1>;
    template <typename T> using const_matrix = const_marray<T, 2>;
    template <typename T> using matrix = marray<T, 2>;
}

#endif
