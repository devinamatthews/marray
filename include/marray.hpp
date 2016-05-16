#ifndef _MARRAY_MARRAY_HPP_
#define _MARRAY_MARRAY_HPP_

#include <type_traits>
#include <array>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <utility>

#include "miterator.hpp"
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

        template<class U>
        struct rebind { typedef aligned_allocator<U, N> other; };
    };

    template <typename T, size_t N, typename U, size_t M>
    bool operator==(const aligned_allocator<T, N>&, const aligned_allocator<U, M>&) { return true; }

    template <typename T, size_t N, typename U, size_t M>
    bool operator!=(const aligned_allocator<T, N>&, const aligned_allocator<U, M>&) { return false; }

    namespace slice
    {
        /*
         * The type all_t specifies a range [0,len_i) for an array
         * dimension i of length len_i (i.e. it selects all of the data along
         * that dimension).
         */
        struct all_t {};
        constexpr all_t all;
    }

    /*
     * The special value uninitialized is used to construct an array which
     * does not default- or value-initialize its elements (useful for avoiding
     * redundant memory operations for scalar types).
     */
    struct uninitialized_t {};
    constexpr uninitialized_t uninitialized;

    struct transpose_t {};
    namespace transpose { constexpr transpose_t T; }

    /*
     * Specifies the layout of the array data.
     */
    enum Layout : int {COLUMN_MAJOR, ROW_MAJOR, DEFAULT=ROW_MAJOR};

    namespace detail
    {
        template <bool IsConst, bool IsView, typename T, unsigned ndim, typename Allocator=std::allocator<T>>
        class marray_base;

        template <bool IsConst, typename T, unsigned ndim, unsigned dim>
        class marray_ref;

        template <bool IsConst, typename T, unsigned ndim, unsigned dim, unsigned newdim>
        class marray_slice;

        /*
         * Represents a part of an array, where the first dim-1 out of ndim
         * dimensions have been indexed into. This type may be implicity converted
         * to an array view or further indexed. This particular reference type
         * is explicitly const, and may only be used to read the original array.
         */
        template <bool IsConst, typename T, unsigned ndim, unsigned dim>
        class marray_ref
        {
            template <bool IsConst_, bool IsView_, typename T_, unsigned ndim_, typename Allocator_>
            friend class marray_base;

            template <bool IsConst_, typename T_, unsigned ndim_, unsigned dim_>
            friend class marray_ref;

            template <bool IsConst_, typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_>
            friend class marray_slice;

            protected:
                typedef marray_base<false, false, T, ndim> base;

                typedef typename base::stride_type stride_type;
                typedef typename base::idx_type idx_type;
                typedef typename base::value_type value_type;
                typedef typename base::pointer pointer;
                typedef typename base::const_pointer const_pointer;
                typedef typename base::reference reference;
                typedef typename base::const_reference const_reference;

                base& array;
                stride_type idx;

                marray_ref(const marray_ref& other) = default;

                template <bool IsConst_, bool IsView_, bool IsConst__=IsConst,
                          typename=typename std::enable_if<IsConst__>::type>
                marray_ref(const marray_base<IsConst_, IsView_, T, ndim>& array, stride_type idx, idx_type i)
                : array(reinterpret_cast<base&>(array)), idx(idx+i*array.stride_[dim-2]) {}

                template <bool IsView_, bool IsConst__=IsConst,
                          typename=typename std::enable_if<!IsConst__>::type>
                marray_ref(marray_base<false, IsView_, T, ndim>& array, stride_type idx, idx_type i)
                : array(reinterpret_cast<base&>(array)), idx(idx+i*array.stride_[dim-2]) {}

                template <bool IsConst__=IsConst,
                          typename=typename std::enable_if<!IsConst__>::type>
                marray_ref(marray_base<false, true, T, ndim>&& array, stride_type idx, idx_type i)
                : array(reinterpret_cast<base&>(array)), idx(idx+i*array.stride_[dim-2]) {}

            public:
                marray_ref& operator=(const marray_ref&) = delete;

                template <unsigned ndim_, bool IsConst_, bool IsConst__=IsConst>
                typename std::enable_if<!IsConst__, marray_ref&>::type
                operator=(const marray_ref<IsConst_, T, ndim_, ndim_-ndim+dim>& other)
                {
                    copy(other, *this);
                    return *this
                }

                template <unsigned ndim_, unsigned newdim_, bool IsConst_, bool IsConst__=IsConst>
                typename std::enable_if<!IsConst__, marray_ref&>::type
                operator=(const marray_slice<IsConst_, T, ndim_, ndim_-ndim+dim+newdim_, newdim_>& other)
                {
                    copy(other, *this);
                    return *this;
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff==0, const T&>::type
                operator[](idx_type i) const
                {
                    return data()[i];
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff==0, T&>::type
                operator[](idx_type i)
                {
                    return data()[i];
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff!=0, marray_ref<true, T, ndim, dim+1>>::type
                operator[](idx_type i) const
                {
                    assert(i >= 0 && i < array.len_[dim-1]);
                    return {array, idx, i};
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff!=0, marray_ref<false, T, ndim, dim+1>>::type
                operator[](idx_type i)
                {
                    assert(i >= 0 && i < array.len_[dim-1]);
                    return {array, idx, i};
                }

                template <typename I, int diff=ndim-dim>
                typename std::enable_if<diff==0, marray_base<true, true, T, 1>>::type
                operator[](const range_t<I>& x) const
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[ndim-1]);
                    return {x.size(), data()+array.stride_[ndim-1]*x.front(), array.stride_[ndim-1]};
                }

                template <typename I, bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff==0, marray_base<false, true, T, 1>>::type
                operator[](const range_t<I>& x)
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[ndim-1]);
                    return {x.size(), data()+array.stride_[ndim-1]*x.front(), array.stride_[ndim-1]};
                }

                template <typename I, int diff=ndim-dim>
                typename std::enable_if<diff!=0, marray_slice<true, T, ndim, dim+1, 1>>::type
                operator[](const range_t<I>& x) const
                {
                    return {array, idx, {}, {}, x};
                }

                template <typename I, bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff!=0, marray_slice<false, T, ndim, dim+1, 1>>::type
                operator[](const range_t<I>& x)
                {
                    return {array, idx, {}, {}, x};
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff==0, marray_base<true, true, T, 1>>::type
                operator[](const slice::all_t& x) const
                {
                    return *this;
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff==0, marray_base<false, true, T, 1>>::type
                operator[](const slice::all_t& x)
                {
                    return *this;
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff!=0, marray_slice<true, T, ndim, dim+1, 1>>::type
                operator[](const slice::all_t& x) const
                {
                    return {array, idx, {}, {}, range(idx_type(), array.len_[dim-1])};
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff!=0, marray_slice<false, T, ndim, dim+1, 1>>::type
                operator[](const slice::all_t& x)
                {
                    return {array, idx, {}, {}, range(idx_type(), array.len_[dim-1])};
                }

                const_pointer data() const
                {
                    return array.data_+idx;
                }

                template <bool IsConst_=IsConst>
                typename std::enable_if<!IsConst_, pointer>::type
                data()
                {
                    return array.data_+idx;
                }

                operator marray_base<true, true, T, ndim-dim+1>() const
                {
                    std::array<idx_type, ndim-dim+1> len;
                    std::array<stride_type, ndim-dim+1> stride;
                    std::copy_n(array.len_.begin(), ndim-dim+1, len.begin());
                    std::copy_n(array.stride_.begin(), ndim-dim+1, stride.begin());
                    return {len, data(), stride};
                }

                template <bool IsConst_=IsConst, typename=typename std::enable_if<!IsConst_>::type>
                operator marray_base<false, true, T, ndim-dim+1>()
                {
                    std::array<idx_type, ndim-dim+1> len;
                    std::array<stride_type, ndim-dim+1> stride;
                    std::copy_n(array.len_.begin(), ndim-dim+1, len.begin());
                    std::copy_n(array.stride_.begin(), ndim-dim+1, stride.begin());
                    return {len, data(), stride};
                }
        };

        template <bool IsConst, typename T, unsigned ndim, unsigned dim>
        marray_base<IsConst, true, T, ndim-dim+1> view(marray_ref<IsConst, T, ndim, dim>& x)
        {
            return x;
        }

        template <bool IsConst, typename T, unsigned ndim, unsigned dim>
        marray_base<IsConst, true, T, ndim-dim+1> view(marray_ref<IsConst, T, ndim, dim>&& x)
        {
            return x;
        }

        /*
         * Represents a part of an array, where the first dim-1 out of ndim
         * dimensions have either been indexed into (i.e. a single value
         * specified for that index) or sliced (i.e. a range of values specified).
         * The parameter newdim specifies how many indices were sliced. The
         * reference may be converted into an array view (of dimension
         * ndim-dim+1+newdim) or further indexed, but may not be used to modify
         * data.
         */
        template <bool IsConst, typename T, unsigned ndim, unsigned dim, unsigned newdim>
        class marray_slice
        {
            template <bool IsConst_, bool IsView_, typename T_, unsigned ndim_, typename Allocator_>
            friend class marray_base;

            template <bool IsConst_, typename T_, unsigned ndim_, unsigned dim_>
            friend class marray_ref;

            template <bool IsConst_, typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_>
            friend class marray_slice;

            protected:
                typedef marray_base<false, false, T, ndim> base;

                typedef typename base::stride_type stride_type;
                typedef typename base::idx_type idx_type;
                typedef typename base::value_type value_type;
                typedef typename base::pointer pointer;
                typedef typename base::const_pointer const_pointer;
                typedef typename base::reference reference;
                typedef typename base::const_reference const_reference;

                base& array;
                stride_type idx;
                std::array<unsigned, newdim> dims;
                std::array<idx_type, newdim> lens;

                marray_slice(const marray_slice& other) = default;

                template <bool IsConst_, bool IsView_, bool IsConst__=IsConst,
                          typename=typename std::enable_if<IsConst__>::type>
                marray_slice(const marray_base<IsConst_, IsView_, T, ndim>& array, stride_type idx,
                             const std::array<unsigned,newdim>& dims,
                             const std::array<idx_type,newdim>& lens, idx_type i)
                : array(reinterpret_cast<base&>(array)), idx(idx+i*array.stride[dim-2]), dims(dims), lens(lens) {}

                template <bool IsView_, bool IsConst__=IsConst,
                          typename=typename std::enable_if<!IsConst__>::type>
                marray_slice(marray_base<false, IsView_, T, ndim>& array, stride_type idx,
                             const std::array<unsigned,newdim>& dims,
                             const std::array<idx_type,newdim>& lens, idx_type i)
                : array(reinterpret_cast<base&>(array)), idx(idx+i*array.stride[dim-2]), dims(dims), lens(lens) {}

                template <bool IsConst__=IsConst,
                          typename=typename std::enable_if<!IsConst__>::type>
                marray_slice(marray_base<false, true, T, ndim>&& array, stride_type idx,
                             const std::array<unsigned,newdim>& dims,
                             const std::array<idx_type,newdim>& lens, idx_type i)
                : array(reinterpret_cast<base&>(array)), idx(idx+i*array.stride[dim-2]), dims(dims), lens(lens) {}

                template <typename I, bool IsConst_, bool IsView_, bool IsConst__=IsConst,
                          typename=typename std::enable_if<IsConst__>::type>
                marray_slice(const marray_base<IsConst_, IsView_, T, ndim>& array, stride_type idx,
                             const std::array<unsigned,newdim-1>& dims,
                             const std::array<idx_type,newdim-1>& lens, const range_t<I>& range_)
                : array(reinterpret_cast<base&>(array)), idx(idx+array.stride_[dim-2]*range_.front())
                {
                    std::copy(dims.begin(), dims.end(), this->dims.begin());
                    this->dims.back() = dim-2;
                    std::copy(lens.begin(), lens.end(), this->lens.begin());
                    this->lens.back() = range_.size();
                }

                template <typename I, bool IsView_, bool IsConst__=IsConst,
                          typename=typename std::enable_if<!IsConst__>::type>
                marray_slice(marray_base<false, IsView_, T, ndim>& array, stride_type idx,
                             const std::array<unsigned,newdim-1>& dims,
                             const std::array<idx_type,newdim-1>& lens, const range_t<I>& range_)
                : array(reinterpret_cast<base&>(array)), idx(idx+array.stride_[dim-2]*range_.front())
                {
                    std::copy(dims.begin(), dims.end(), this->dims.begin());
                    this->dims.back() = dim-2;
                    std::copy(lens.begin(), lens.end(), this->lens.begin());
                    this->lens.back() = range_.size();
                }

                template <typename I, bool IsConst__=IsConst,
                          typename=typename std::enable_if<!IsConst__>::type>
                marray_slice(marray_base<false, true, T, ndim>&& array, stride_type idx,
                             const std::array<unsigned,newdim-1>& dims,
                             const std::array<idx_type,newdim-1>& lens, const range_t<I>& range_)
                : array(reinterpret_cast<base&>(array)), idx(idx+array.stride_[dim-2]*range_.front())
                {
                    std::copy(dims.begin(), dims.end(), this->dims.begin());
                    this->dims.back() = dim-2;
                    std::copy(lens.begin(), lens.end(), this->lens.begin());
                    this->lens.back() = range_.size();
                }

            public:
                marray_slice& operator=(const marray_slice&) = delete;

                template <unsigned ndim_, bool IsConst_, bool IsConst__=IsConst>
                typename std::enable_if<!IsConst__, marray_slice&>::type
                operator=(const marray_ref<IsConst_, T, ndim_, ndim_-ndim+dim-newdim>& other)
                {
                    copy(other, *this);
                    return *this;
                }

                template <unsigned ndim_, unsigned newdim_, bool IsConst_, bool IsConst__=IsConst>
                typename std::enable_if<!IsConst__, marray_slice&>::type
                operator=(const marray_slice<IsConst_, T, ndim_, ndim_-ndim+dim-newdim+newdim_, newdim_>& other)
                {
                    copy(other, *this);
                    return *this;
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff==0, marray_base<true, true, T, newdim>>::type
                operator[](idx_type i) const
                {
                    assert(i >= 0 && i < array.len_[ndim-1]);
                    std::array<stride_type, newdim> strides;
                    for (unsigned i = 0;i < newdim;i++)
                    {
                        strides[i] = array.stride_[dims[i]];
                    }
                    return {lens, data()+i*array.stride_[ndim-1], strides};
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff==0, marray_base<false, true, T, newdim>>::type
                operator[](idx_type i)
                {
                    assert(i >= 0 && i < array.len_[ndim-1]);
                    std::array<stride_type, newdim> strides;
                    for (unsigned i = 0;i < newdim;i++)
                    {
                        strides[i] = array.stride_[dims[i]];
                    }
                    return {lens, data()+i*array.stride_[ndim-1], strides};
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff!=0, marray_slice<true, T, ndim, dim+1, newdim>>::type
                operator[](idx_type i) const
                {
                    assert(i >= 0 && i < array.len[dim-1]);
                    return {array, idx, dims, lens, i};
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff!=0, marray_slice<false, T, ndim, dim+1, newdim>>::type
                operator[](idx_type i)
                {
                    assert(i >= 0 && i < array.len[dim-1]);
                    return {array, idx, dims, lens, i};
                }

                template <typename I, int diff=ndim-dim>
                typename std::enable_if<diff==0, marray_base<true, true, T, newdim+1>>::type
                operator[](const range_t<I>& x) const
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[ndim-1]);
                    std::array<idx_type, newdim+1> newlens;
                    std::array<stride_type, newdim+1> strides;
                    for (unsigned i = 0;i < newdim;i++)
                    {
                        newlens[i] = lens[i];
                        strides[i] = array.stride_[dims[i]];
                    }
                    newlens[newdim] = x.size();
                    strides[newdim] = array.stride_[ndim-1];
                    return {newlens, data()+array.stride_[ndim-1]*x.front(), strides};
                }

                template <typename I, bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff==0, marray_base<false, true, T, newdim+1>>::type
                operator[](const range_t<I>& x)
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[ndim-1]);
                    std::array<idx_type, newdim+1> newlens;
                    std::array<stride_type, newdim+1> strides;
                    for (unsigned i = 0;i < newdim;i++)
                    {
                        newlens[i] = lens[i];
                        strides[i] = array.stride_[dims[i]];
                    }
                    newlens[newdim] = x.size();
                    strides[newdim] = array.stride_[ndim-1];
                    return {newlens, data()+array.stride_[ndim-1]*x.front(), strides};
                }

                template <typename I, int diff=ndim-dim>
                typename std::enable_if<diff!=0, marray_slice<true, T, ndim, dim+1, newdim+1>>::type
                operator[](const range_t<I>& x) const
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[dim-1]);
                    return {array, idx, dims, lens, x};
                }

                template <typename I, bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff!=0, marray_slice<false, T, ndim, dim+1, newdim+1>>::type
                operator[](const range_t<I>& x)
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[dim-1]);
                    return {array, idx, dims, lens, x};
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff==0, marray_base<true, true, T, newdim+1>>::type
                operator[](const slice::all_t& x) const
                {
                    return *this;
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff==0, marray_base<false, true, T, newdim+1>>::type
                operator[](const slice::all_t& x)
                {
                    return *this;
                }

                template <int diff=ndim-dim>
                typename std::enable_if<diff!=0, marray_slice<true, T, ndim, dim+1, newdim+1>>::type
                operator[](const slice::all_t& x) const
                {
                    return {array, idx, dims, lens, range(idx_type(), array.len_[dim-1])};
                }

                template <bool IsConst_=IsConst, int diff=ndim-dim>
                typename std::enable_if<!IsConst && diff!=0, marray_slice<false, T, ndim, dim+1, newdim+1>>::type
                operator[](const slice::all_t& x)
                {
                    return {array, idx, dims, lens, range(idx_type(), array.len_[dim-1])};
                }

                const_pointer data() const
                {
                    return array.data_+idx;
                }

                template <bool IsConst_=IsConst>
                typename std::enable_if<!IsConst, pointer>::type
                data()
                {
                    return array.data_+idx;
                }

                operator marray_base<true, true, T, ndim+newdim-dim+1>() const
                {
                    std::array<idx_type, ndim+newdim-dim+1> len;
                    std::array<stride_type, ndim+newdim-dim+1> stride;
                    for (unsigned i = 0;i < newdim;i++)
                    {
                        len[i] = lens[i];
                        stride[i] = array.stride_[dims[i]];
                    }
                    std::copy_n(array.len_.begin()+dim-1, ndim-dim+1, len.begin()+newdim);
                    std::copy_n(array.stride_.begin()+dim-1, ndim-dim+1, stride.begin()+newdim);
                    return {len, data(), stride};
                }

                template <bool IsConst_=IsConst, typename=typename std::enable_if<!IsConst>::type>
                operator marray_base<false, true, T, ndim+newdim-dim+1>()
                {
                    std::array<idx_type, ndim+newdim-dim+1> len;
                    std::array<stride_type, ndim+newdim-dim+1> stride;
                    for (unsigned i = 0;i < newdim;i++)
                    {
                        len[i] = lens[i];
                        stride[i] = array.stride_[dims[i]];
                    }
                    std::copy_n(array.len_.begin()+dim-1, ndim-dim+1, len.begin()+newdim);
                    std::copy_n(array.stride_.begin()+dim-1, ndim-dim+1, stride.begin()+newdim);
                    return {len, data(), stride};
                }
        };

        template <bool IsConst, typename T, unsigned ndim, unsigned dim, unsigned newdim>
        marray_base<IsConst, true, T, ndim+newdim-dim+1> view(marray_slice<IsConst, T, ndim, dim, newdim>& x)
        {
            return x;
        }

        template <bool IsConst, typename T, unsigned ndim, unsigned dim, unsigned newdim>
        marray_base<IsConst, true, T, ndim+newdim-dim+1> view(marray_slice<IsConst, T, ndim, dim, newdim>&& x)
        {
            return x;
        }

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
        struct marray_type;

        template <typename T, unsigned N>
        struct marray_type<false, T, N> { typedef marray_base<false, true, T, N> type; };
        
        template <typename T, unsigned N>
        struct marray_type<true, T, N> { typedef marray_base<true, true, T, N> type; };

        template <typename T>
        struct marray_type<false, T, 0> { typedef T& type; };

        template <typename T>
        struct marray_type<true, T, 0> { typedef const T& type; };

        /*
         * These helper structs determine the number of sliced dimensions based
         * on a set of indices which may be either an integral type (no slicing)
         * or a range (slicing).
         */
        template <typename Index>
        struct is_slice : std::false_type {};

        template <typename I>
        struct is_slice<range_t<I>> : std::true_type {};

        template <typename... Indices>
        struct num_slices;

        template <>
        struct num_slices<>
        {
            static constexpr unsigned N = 0;
        };

        template <typename Index, typename... Indices>
        struct num_slices<Index, Indices...>
        {
            static constexpr unsigned N = num_slices<Indices...>::N + is_slice<Index>::value;
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
            RT operator()(ArrayOrSlice&& x, Index idx0, Indices... idx)
            {
                return get_slice<RT, Indices...>()(x[idx0], idx...);
            }
        };

        template <typename RT, typename Index>
        struct get_slice<RT, Index>
        {
            template <typename ArrayOrSlice>
            RT operator()(ArrayOrSlice&& x, Index idx0)
            {
                return x[idx0];
            }
        };

        template <unsigned N, unsigned I, template <typename...> class Condition, typename... Args>
        struct are_integral_helper : std::false_type {};

        template <unsigned N, unsigned I, template <typename...> class Condition, typename Arg, typename... Args>
        struct are_integral_helper<N, I, Condition, Arg, Args...>
        : std::conditional<std::is_integral<Arg>::value,
                           are_integral_helper<N, I+1, Condition, Args...>,
                           std::false_type>::type {};

        template <unsigned N, template <typename...> class Condition, typename Arg, typename... Args>
        struct are_integral_helper<N, N, Condition, Arg, Args...>
        : std::conditional<std::is_integral<Arg>::value,
                           Condition<Args...>,
                           std::false_type>::type {};

        template <unsigned N, template <typename...> class Condition, typename... Args>
        using are_integral =
            typename std::conditional<N,
                                      are_integral_helper<N, 1, Condition, Args...>,
                                      Condition<Args...>>::type;

        template <typename... Args>
        struct are_empty : std::integral_constant<bool, (sizeof...(Args) == 0)> {};

        /*
         * Helper classes to set the len_ and stride_ arrays of an marray.
         */
        template <unsigned N, unsigned I, typename Tail, typename Arg, typename... Args>
        struct set_len_helper
        {
            template <typename Array>
            set_len_helper(Array&& array, std::array<size_t, N>& len, const Arg& arg, const Args&... args)
            {
                len[I-1] = arg;
                set_len_helper<N, I+1, Tail, Args...> tmp(array, len, args...);
            }
        };

        template <unsigned N, typename Tail, typename Arg, typename... Args>
        struct set_len_helper<N, N, Tail, Arg, Args...>
        {
            template <typename Array>
            set_len_helper(Array&& array, std::array<size_t, N>& len, const Arg& arg, const Args&... args)
            {
                len[N-1] = arg;
                Tail tmp(array, len, args...);
            }
        };

        template <unsigned N, typename Tail, typename... Args>
        using set_len = set_len_helper<N, 1, Tail, Args...>;

        /*
         * Helper class to determine if the arguments Args... are of one of the
         * forms which are acceptable for passing to marray<T, ndim>::reset.
         */
        template <typename T, unsigned ndim>
        struct are_reset_args
        {
            template <typename... Args> struct direct_helper;

            template <typename... Args>
            using direct = std::integral_constant<bool,
                are_integral<ndim, direct_helper,
                             typename std::decay<Args>::type...>::value>;

            template <typename... Args> struct const_view_helper;

            template <typename... Args>
            using const_view = std::integral_constant<bool,
                are_integral<ndim, const_view_helper,
                             typename std::decay<Args>::type...>::value>;

            template <typename... Args> struct view_helper;

            template <typename... Args>
            using view = std::integral_constant<bool,
                are_integral<ndim, view_helper,
                             typename std::decay<Args>::type...>::value>;
        };

        template <typename T, unsigned ndim>
        template <typename... Args>
        struct are_reset_args<T,ndim>::direct_helper : std::integral_constant<bool,!sizeof...(Args)> {};

        template <typename T, unsigned ndim>
        template <typename Arg>
        struct are_reset_args<T,ndim>::direct_helper<Arg>
        : std::integral_constant<bool, std::is_convertible<Arg,T>::value ||
                                       std::is_same<Arg,uninitialized_t>::value ||
                                       std::is_same<Arg,Layout>::value> {};

        template <typename T, unsigned ndim>
        template <typename Arg1, typename Arg2>
        struct are_reset_args<T,ndim>::direct_helper<Arg1, Arg2>
        : std::integral_constant<bool, (std::is_convertible<Arg1,T>::value ||
                                        std::is_same<Arg1,uninitialized_t>::value) &&
                                       std::is_same<Arg2,Layout>::value> {};

        template <typename T, unsigned ndim>
        template <typename... Args>
        struct are_reset_args<T,ndim>::const_view_helper : std::false_type {};

        template <typename T, unsigned ndim>
        template <typename Arg>
        struct are_reset_args<T,ndim>::const_view_helper<Arg> : std::is_convertible<Arg,const T*> {};

        template <typename T, unsigned ndim>
        template <typename Arg>
        struct are_reset_args<T,ndim>::const_view_helper<Arg, Layout> : std::is_convertible<Arg,const T*> {};

        template <typename T, unsigned ndim>
        template <typename Arg, typename... Args>
        struct are_reset_args<T,ndim>::const_view_helper<Arg, Args...>
        : std::integral_constant<bool, std::is_convertible<Arg,const T*>::value &&
                                       are_integral<ndim, are_empty, Args...>::value> {};

        template <typename T, unsigned ndim>
        template <typename... Args>
        struct are_reset_args<T,ndim>::view_helper : std::false_type {};

        template <typename T, unsigned ndim>
        template <typename Arg>
        struct are_reset_args<T,ndim>::view_helper<Arg> : std::is_convertible<Arg,T*> {};

        template <typename T, unsigned ndim>
        template <typename Arg>
        struct are_reset_args<T,ndim>::view_helper<Arg, Layout> : std::is_convertible<Arg,T*> {};

        template <typename T, unsigned ndim>
        template <typename Arg, typename... Args>
        struct are_reset_args<T,ndim>::view_helper<Arg, Args...>
        : std::integral_constant<bool, std::is_convertible<Arg,T*>::value &&
                                       are_integral<ndim, are_empty, Args...>::value> {};

        struct reset_helper
        {
            template <typename T, unsigned ndim>
            reset_helper(marray<T, ndim>& array, std::array<size_t, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
            {
                array.reset(len, val, layout);
            }

            template <typename T, unsigned ndim>
            reset_helper(marray<T, ndim>& array, std::array<size_t, ndim>& len, const uninitialized_t& u, Layout layout = DEFAULT)
            {
                array.reset(len, u, layout);
            }

            template <typename T, unsigned ndim>
            reset_helper(const_marray_view<T, ndim>& array, std::array<size_t, ndim>& len, const T* ptr, Layout layout = DEFAULT)
            {
                array.reset(len, const_cast<T*>(ptr), layout);
            }

            template <typename T, unsigned ndim, typename... Args>
            reset_helper(const_marray_view<T, ndim>& array, std::array<size_t, ndim>& len, const T* ptr, Args... args)
            {
                array.reset(len, const_cast<T*>(ptr), make_array(args...));
            }
        };

        template <typename T, unsigned ndim, typename... Args>
        void reset(const_marray_view<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<ndim, reset_helper, Args...> tmp(array, len, args...);
        }

        template <typename T, unsigned ndim, typename... Args>
        void reset(marray_view<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<ndim, reset_helper, Args...> tmp(array, len, args...);
        }

        template <typename T, unsigned ndim, typename... Args>
        void reset(marray<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<ndim, reset_helper, Args...> tmp(array, len, args...);
        }

        /*
         * Helper class to determine if the arguments Args... are of one of the
         * forms which are acceptable for passing to marray<T, ndim>::resize.
         */
        template <typename T, unsigned ndim>
        struct are_resize_args
        {
            template <typename... Args> struct direct_helper;

            template <typename... Args>
            using direct =  std::integral_constant<bool,
                 are_integral<ndim, direct_helper,
                              typename std::decay<Args>::type...>::value>;
        };

        template <typename T, unsigned ndim>
        template <typename... Args>
        struct are_resize_args<T,ndim>::direct_helper : std::integral_constant<bool,!sizeof...(Args)> {};

        template <typename T, unsigned ndim>
        template <typename Arg>
        struct are_resize_args<T,ndim>::direct_helper<Arg>
        : std::integral_constant<bool, std::is_convertible<Arg,T>::value> {};

        struct resize_helper
        {
            template <typename T, unsigned ndim>
            resize_helper(marray<T, ndim>& array, std::array<size_t, ndim>& len, const T& val=T())
            {
                array.resize(len, val);
            }
        };

        template <typename T, unsigned ndim, typename... Args>
        void resize(marray<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<ndim, resize_helper, Args...> tmp(array, len, args...);
        }

        template <bool IsConst, bool IsView, typename T, unsigned ndim, typename Allocator>
        class marray_base
        {
            static_assert(ndim > 0, "0-dimensional marrays are not allowed.");

            template <bool IsConst_, bool IsView_, typename T_, unsigned ndim_, typename Allocator_>
            friend class marray_base;

            template <bool IsConst_, typename T_, unsigned ndim_, unsigned dim_>
            friend class marray_ref;

            template <bool IsConst_, typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_>
            friend class marray_slice;

            public:
                typedef unsigned idx_type;
                typedef size_t size_type;
                typedef ptrdiff_t stride_type;
                typedef T value_type;
                typedef T* pointer;
                typedef const T* const_pointer;
                typedef T& reference;
                typedef const T& const_reference;

            protected:
                struct mem_t : Allocator
                {
                    typedef std::allocator_traits<Allocator> alloc_traits;

                    mem_t(const mem_t&)  = delete;

                    mem_t& operator=(const mem_t&) = delete;

                    mem_t(Allocator alloc) : Allocator(std::move(alloc)) {}

                    ~mem_t() { clear(); }

                    void allocate(size_type size, uninitialized_t u)
                    {
                        clear();
                        size_ = size;
                        is_alloced_ = true;
                        data_ = alloc_traits::allocate(*this, size);
                    }

                    void allocate(size_type size, const T& value=T())
                    {
                        allocate(size, uninitialized);
                        for (size_type i = 0;i < size;i++)
                        {
                            alloc_traits::construct(*this, data_+i, value);
                        }
                    }

                    void set(pointer data)
                    {
                        clear();
                        data_ = data;
                    }

                    void set(mem_t&& other)
                    {
                        using std::swap;
                        clear();
                        swap(other.is_alloced_, is_alloced_);
                        swap(other.data_, data_);
                        swap(other.size_, size_);
                    }

                    void set(const mem_t& other)
                    {
                        clear();
                        allocate(other.size_, uninitialized);
                        std::uninitialized_copy_n(other.data_, size_, data_);
                    }

                    void clear()
                    {
                        if (is_alloced_ && data_)
                        {
                            for (size_type i = 0;i < size_;i++)
                            {
                                alloc_traits::destroy(*this, data_+i);
                            }
                            alloc_traits::deallocate(*this, data_, size_);
                        }
                        data_ = nullptr;
                        is_alloced_ = false;
                    }

                    pointer data_ = nullptr;
                    size_t size_ = 0;
                    bool is_alloced_ = false;
                } mem_;

                std::array<idx_type,ndim> len_ = {};
                std::array<stride_type,ndim> stride_ = {};
                Layout layout_ = DEFAULT;

                marray_base(const marray_base& other) = delete;

                marray_base(marray_base&& other) = delete;

                marray_base& operator=(const marray_base& other) = delete;

                marray_base& operator=(marray_base&& other) = delete;

            public:
                template <typename U>
                static std::array<stride_type, ndim> default_strides(const std::array<U, ndim>& len, Layout layout=DEFAULT)
                {
                    std::array<stride_type, ndim> stride;

                    if (layout == ROW_MAJOR)
                    {
                        stride[ndim-1] = 1;
                        for (unsigned i = ndim-1;i > 0;i--)
                        {
                            stride[i-1] = stride[i]*len[i];
                        }
                    }
                    else
                    {
                        stride[0] = 1;
                        for (unsigned i = 1;i < ndim;i++)
                        {
                            stride[i] = stride[i-1]*len[i-1];
                        }
                    }

                    return stride;
                }

                template <typename... Args>
                static typename std::enable_if<detail::are_integral<ndim, detail::are_empty, typename std::decay<Args>::type...>::value>::type
                default_strides(Args... args)
                {
                    return default_strides(make_array<idx_type>(args...));
                }

                marray_base() {}

                template <bool IsConst_, bool IsView_, bool IsConst__=IsConst, bool IsView__=IsView,
                          typename=typename std::enable_if<IsConst__ && IsView__>::type>
                marray_base(const marray_base<IsConst_, IsView_, T, ndim>& other)
                {
                    view(other);
                }

                template <bool IsView_, bool IsConst__=IsConst, bool IsView__=IsView,
                          typename=typename std::enable_if<!IsConst__ && IsView__>::type>
                marray_base(marray_base<false, IsView_, T, ndim>& other)
                {
                    view(other);
                }

                template <bool IsConst__=IsConst, bool IsView__=IsView,
                          typename=typename std::enable_if<!IsConst__ && IsView__>::type>
                marray_base(marray_base<false, true, T, ndim>&& other)
                {
                    view(other);
                }

                template <typename U, bool IsConst__=IsConst, bool IsView__=IsView,
                          typename=typename std::enable_if<IsConst__ && IsView__>::type>
                marray_base(const std::array<U, ndim>& len, const_pointer ptr, Layout layout = DEFAULT)
                {
                    reset(len, ptr, layout);
                }

                template <typename U, typename V, bool IsConst__=IsConst, bool IsView__=IsView,
                          typename=typename std::enable_if<IsConst__ && IsView__>::type>
                marray_base(const std::array<U, ndim>& len, const_pointer ptr, const std::array<V, ndim>& stride)
                {
                    reset(len, ptr, stride);
                }

                template <typename U, bool IsView__=IsView,
                          typename=typename std::enable_if<IsView__>::type>
                marray_base(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
                {
                    reset(len, ptr, layout);
                }

                template <typename U, typename V, bool IsView__=IsView,
                          typename=typename std::enable_if<IsView__>::type>
                marray_base(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
                {
                    reset(len, ptr, stride);
                }

                template <bool IsConst_, bool IsView_, bool IsView__=IsView,
                          typename=typename std::enable_if<!IsView__>::type>
                marray_base(const marray_base<IsConst_, IsView_, T, ndim>& other, Layout layout=DEFAULT)
                : layout_(layout)
                {
                    copy(other, layout);
                }

                template <bool IsConst_, bool IsView__=IsView,
                          typename=typename std::enable_if<!IsView__>::type>
                marray_base(marray_base<IsConst_, false, T, ndim>&& other, Layout layout=DEFAULT)
                : layout_(layout)
                {
                    copy(std::move(other), layout);
                }

                template <typename U>
                marray_base(const std::array<U, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
                : layout_(layout)
                {
                    reset(len, val, layout);
                }

                template <typename U>
                marray_base(const std::array<U, ndim>& len, uninitialized_t u, Layout layout = DEFAULT)
                : layout_(layout)
                {
                    reset(len, u, layout);
                }

                template <typename... Args, typename=
                    typename std::enable_if<(IsConst && IsView &&
                                             detail::are_reset_args<T, ndim>
                                                 ::template const_view<Args...>::value) ||
                                            (!IsConst && IsView &&
                                             detail::are_reset_args<T, ndim>
                                                 ::template view<Args...>::value) ||
                                            (!IsView &&
                                             detail::are_reset_args<T, ndim>
                                                 ::template direct<Args...>::value)
                                           >::type>
                explicit marray_base(Args&&... args)
                {
                    reset(std::forward<Args>(args)...);
                }

                template <bool IsConst_, bool IsView_, bool IsConst__=IsConst>
                typename std::enable_if<!IsConst__, marray_base&>::type
                operator=(const marray_base<IsConst_, IsView_, T,ndim>& other)
                {
                    copy(other, *this);
                    return *this;
                }

                template <bool IsConst_, bool IsView_, bool IsConst__=IsConst, bool IsView__=IsView>
                typename std::enable_if<IsConst__ && IsView__>::type
                view(const marray_base<IsConst_, IsView, T, ndim>& other)
                {
                    mem_.set(other.mem_.data_);
                    len_ = other.len_;
                    stride_ = other.stride_;
                }

                template <bool IsView_, bool IsConst__=IsConst, bool IsView__=IsView>
                typename std::enable_if<!IsConst__ && IsView__>::type
                view(marray_base<false, IsView, T, ndim>& other)
                {
                    mem_.set(other.mem_.data_);
                    len_ = other.len_;
                    stride_ = other.stride_;
                }

                template <bool IsConst__=IsConst, bool IsView__=IsView>
                typename std::enable_if<!IsConst__ && IsView__>::type
                view(marray_base<false, true, T, ndim>&& other)
                {
                    mem_.set(other.mem_.data_);
                    len_ = other.len_;
                    stride_ = other.stride_;
                }

                template <bool IsConst_, bool IsView_, bool IsView__=IsView>
                typename std::enable_if<!IsView__>::type
                copy(const marray_base<IsConst_, IsView_, T, ndim>& other, Layout layout=DEFAULT)
                {
                    if (std::is_scalar<T>::value)
                    {
                        reset(other.len_, uninitialized, layout);
                    }
                    else
                    {
                        reset(other.len_, T(), layout);
                    }

                    *this = other;
                }

                template <bool IsView__=IsView>
                typename std::enable_if<!IsView__>::type
                copy(marray_base<false, false, T, ndim>&& other)
                {
                    other.mem_.set(std::move(other.mem_));
                    len_ = other.len_;
                    stride_ = other.stride_;
                    other.len_.fill(0);
                }

                void reset()
                {
                    mem_.clear();
                    len_.fill(0);
                    stride_.fill(0);
                }

                template <typename U, bool IsView_=IsView>
                typename std::enable_if<!IsView_>::type
                reset(const std::array<U, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
                {
                    reset(len, uninitialized, layout);
                    std::uninitialized_fill_n(mem_.data_, mem_.size_, val);
                }

                template <typename U, bool IsView_=IsView>
                typename std::enable_if<!IsView_>::type
                reset(const std::array<U, ndim>& len, uninitialized_t u, Layout layout = DEFAULT)
                {
                    layout_ = layout;
                    std::copy_n(len.begin(), ndim, len_.begin());
                    stride_ = default_strides(len_, layout_);
                    mem_.allocate(std::accumulate(len_.begin(), len_.end(), size_t(1), std::multiplies<size_t>()), u);
                }

                template <typename U, bool IsView_=IsView>
                typename std::enable_if<IsView_>::type
                reset(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
                {
                    reset(len, ptr, default_strides(len, layout));
                }

                template <typename U, bool IsConst_=IsConst, bool IsView_=IsView>
                typename std::enable_if<IsConst_ && IsView_>::type
                reset(const std::array<U, ndim>& len, const_pointer ptr, Layout layout = DEFAULT)
                {
                    reset(len, ptr, default_strides(len, layout));
                }

                template <typename U, typename V, bool IsView_=IsView>
                typename std::enable_if<IsView_>::type
                reset(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
                {
                    mem_.set(ptr);
                    std::copy_n(len.begin(), ndim, len_.begin());
                    std::copy_n(stride.begin(), ndim, stride_.begin());
                }

                template <typename U, typename V, bool IsConst_=IsConst, bool IsView_=IsView>
                typename std::enable_if<IsConst_ && IsView_>::type
                reset(const std::array<U, ndim>& len, const_pointer ptr, const std::array<V, ndim>& stride)
                {
                    mem_.set(ptr);
                    std::copy_n(len.begin(), ndim, len_.begin());
                    std::copy_n(stride.begin(), ndim, stride_.begin());
                }

                template <typename... Args, typename=
                    typename std::enable_if<(IsConst && IsView &&
                                             detail::are_reset_args<T, ndim>
                                                 ::template const_view<Args...>::value) ||
                                            (!IsConst && IsView &&
                                             detail::are_reset_args<T, ndim>
                                                 ::template view<Args...>::value) ||
                                            (!IsView &&
                                             detail::are_reset_args<T, ndim>
                                                 ::template direct<Args...>::value)
                                           >::type>
                void reset(Args&&... args)
                {
                    detail::reset(*this, std::forward<Args>(args)...);
                }

                template <typename U, typename V, bool IsView_=IsView>
                typename std::enable_if<!IsView_>::type
                resize(const std::array<U, ndim>& len, const T& val=T())
                {
                    marray_base a(std::move(*this));
                    reset(len, val, layout_);
                    marray_base<true, true, T, ndim> b(*this);

                    /*
                     * It is OK to change the geometry of 'a' even if it is not
                     * a view since it is about to go out of scope.
                     */
                    for (unsigned i = 0;i < ndim;i++)
                    {
                        a.len_[i] = b.len_[i] = std::min(a.len_[i], b.len_[i]);
                    }

                    copy(a, b);
                }

                template <typename... Args, typename=
                    typename std::enable_if<!IsView &&
                                            detail::are_resize_args<T, ndim>
                                                ::template direct<Args...>::value
                                           >::type>
                void resize(Args&&... args)
                {
                    detail::resize(*this, std::forward<Args>(args)...);
                }

                template <bool IsView_=IsView, unsigned ndim_=ndim>
                typename std::enable_if<!IsView_ && ndim_==1>::type
                push_back(const T& x)
                {
                    resize(len_[0]+1);
                    back() = x;
                }

                template <bool IsView_=IsView, unsigned ndim_=ndim>
                typename std::enable_if<!IsView_ && ndim_==1>::type
                pop_back()
                {
                    resize(len_[0]-1);
                }

                template <bool IsConst_, bool IsView_, bool IsView__=IsView, unsigned ndim_=ndim>
                typename std::enable_if<!IsView__ && ndim_!=1>::type
                push_back(unsigned dim, const marray_base<IsConst_, IsView_, T, ndim-1>& x)
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
                    this->back() = x;
                }

                template <bool IsView__=IsView, unsigned ndim_=ndim>
                typename std::enable_if<!IsView__ && ndim_!=1>::type
                pop_back(unsigned dim)
                {
                    assert(dim < ndim);
                    assert(len_[dim] > 0);

                    std::array<idx_type, ndim> len = len_;
                    len[dim]--;
                    resize(len);
                }

                template <typename U>
                marray_base<IsConst, true, T, ndim> permute(const std::array<U, ndim>& perm)
                {
                    marray_base<IsConst, true, T, ndim> view;

                    view.mem_.set(mem_.data_);

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
                typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, marray_base<IsConst, true, T, ndim>>::type
                permute(Args&&... args)
                {
                    return permute(make_array(std::forward<Args>(args)...));
                }

                template <typename U>
                marray_base<IsConst, true, T, ndim> permute(const std::array<U, ndim>& perm)
                {
                    marray_base<IsConst, true, T, ndim> view;

                    view.mem_.set(mem_.data_);

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
                typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, marray_base<IsConst, true, T, ndim>>::type
                permute(Args&&... args)
                {
                    return permute(make_array(std::forward<Args>(args)...));
                }

                template <bool IsView_=IsView>
                typename std::enable_if<!IsView_>::type
                rotate_dim(unsigned dim, idx_type shift)
                {
                    assert(dim < ndim);

                    idx_type n = len_[dim];
                    stride_type s = stride_[dim];

                    shift = shift%n;
                    if (shift < 0) shift += n;

                    if (shift == 0) return;

                    std::array<idx_type, ndim-1> sublen;
                    std::array<stride_type, ndim-1> substride;

                    std::copy_n(len_.begin(), dim, sublen.begin());
                    std::copy_n(len_.begin()+dim+1, ndim-dim-1, sublen.begin()+dim);

                    std::copy_n(stride_.begin(), dim, substride.begin());
                    std::copy_n(stride_.begin()+dim+1, ndim-dim-1, substride.begin()+dim);

                    pointer p = mem_.data_;
                    miterator<idx_type, stride_type, ndim> it(sublen, substride);
                    while (it.next(p))
                    {
                        pointer a = p;
                        pointer b = p+(shift-1)*s;
                        for (idx_type i = 0;i < shift/2;i++)
                        {
                            std::swap(*a, *b);
                            a += s;
                            b -= s;
                        }

                        a = p+shift*s;
                        b = p+(n-1)*s;
                        for (idx_type i = 0;i < (n-shift)/2;i++)
                        {
                            std::swap(*a, *b);
                            a += s;
                            b -= s;
                        }

                        a = p;
                        b = p+(n-1)*s;
                        for (idx_type i = 0;i < n/2;i++)
                        {
                            std::swap(*a, *b);
                            a += s;
                            b -= s;
                        }
                    }
                }

                template <typename U, bool IsView_=IsView>
                typename std::enable_if<!IsView_>::type
                rotate(const std::array<U, ndim>& shift)
                {
                    for (unsigned dim = 0;dim < ndim;dim++)
                    {
                        rotate_dim(dim, shift[dim]);
                    }
                }

                template <typename... Args>
                typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value && !IsView>::type
                rotate(Args&&... args)
                {
                    rotate(make_array(std::forward<Args>(args)...));
                }

                template <unsigned ndim_=ndim>
                typename std::enable_if<ndim_==1, const_reference>::type
                front() const
                {
                    assert(length() > 0);
                    return mem_.data_[0];
                }

                template <bool IsConst_=IsConst, unsigned ndim_=ndim>
                typename std::enable_if<!IsConst && ndim_==1, reference>::type
                front()
                {
                    assert(length() > 0);
                    return mem_.data_[0];
                }

                template <unsigned ndim_=ndim>
                typename std::enable_if<ndim_==1, const_reference>::type
                front(unsigned dim) const
                {
                    assert(dim == 0);
                    return front();
                }

                template <bool IsConst_=IsConst, unsigned ndim_=ndim>
                typename std::enable_if<!IsConst && ndim_==1, reference>::type
                front(unsigned dim)
                {
                    assert(dim == 0);
                    return front();
                }

                template <unsigned ndim_=ndim>
                typename std::enable_if<ndim_!=1>::type>
                marray_base<T, ndim-1> front(unsigned dim)
                {
                    assert(dim < ndim);
                    assert(len_[dim] > 0);

                    marray_view<T, ndim-1> view;
                    view.mem_.set(mem_.data_);
                    copy_n(stride_.begin(), dim, view.stride_.begin());
                    copy_n(stride_.begin()+dim+1, ndim-dim-1, view.stride_.begin()+dim);
                    copy_n(len_.begin(), dim, view.len_.begin());
                    copy_n(len_.begin()+dim+1, ndim-dim-1, view.len_.begin()+dim);

                    return view;
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                reference back()
                {
                    assert(length() > 0);
                    return mem_.data_[(length()-1)*stride()];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                reference back(unsigned dim)
                {
                    assert(dim == 0);
                    return back();
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
                marray_view<T, ndim-1> back(unsigned dim)
                {
                    marray_view<T, ndim-1> view = front(dim);
                    view.data_ += (len_[dim]-1)*stride_[dim];
                    return view;
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                const_reference front() const
                {
                    assert(length() > 0);
                    return mem_.data_[0];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                const_reference front(unsigned dim) const
                {
                    assert(dim == 0);
                    return front();
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
                const_marray_view<T, ndim-1> front(unsigned dim) const
                {
                    return const_cast<marray_base&>(*this).front(dim);
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                const_reference back() const
                {
                    assert(length() > 0);
                    return mem_.data_[(length()-1)*stride()];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                const_reference back(unsigned dim) const
                {
                    assert(dim == 0);
                    return back();
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
                const_marray_view<T, ndim-1> back(unsigned dim) const
                {
                    return const_cast<marray_base&>(*this).back(dim);
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                reference operator[](idx_type i)
                {
                    assert(i < length());
                    return mem_.data_[i*stride()];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
                marray_ref<T, ndim, 2> operator[](idx_type i)
                {
                    assert(i < len_[0]);
                    return const_marray_ref<T, ndim, 2>(*this, (stride_type)0, i);
                }

                template <typename I>
                marray_slice<T, ndim, 2, 1> operator[](const range_t<I>& x)
                {
                    assert(x.front() >= 0 && x.back() <= len_[0]);
                    return const_marray_slice<T, ndim, 2, 1>(*this, (stride_type)0, {}, {}, x);
                }

                marray_slice<T, ndim, 2, 1> operator[](slice::all_t x)
                {
                    return const_marray_slice<T, ndim, 2, 1>(*this, (stride_type)0, {}, {}, range(idx_type(), len_[0]));
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                const_reference operator[](idx_type i) const
                {
                    assert(i < length());
                    return mem_.data_[i*stride()];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
                const_marray_ref<T, ndim, 2> operator[](idx_type i) const
                {
                    return const_cast<marray_base&>(*this)[i];
                }

                template <typename I>
                const_marray_slice<T, ndim, 2, 1> operator[](const range_t<I>& x) const
                {
                    return const_cast<marray_base&>(*this)[x];
                }

                const_marray_slice<T, ndim, 2, 1> operator[](slice::all_t x) const
                {
                    return const_cast<marray_base&>(*this)[x];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                reference operator()(idx_type i)
                {
                    return (*this)[i];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                const_reference operator()(idx_type i) const
                {
                    return (*this)[i];
                }

                template <typename... Indices>
                typename std::enable_if<(sizeof...(Indices) == ndim) && (sizeof...(Indices) > 1),
                                        typename detail::return_type<false, T, Indices...>::type>::type
                operator()(Indices... idx)
                {
                    return detail::get_slice<typename detail::return_type<true, T, Indices...>::type, Indices...>()(*this, idx...);
                }

                template <typename... Indices>
                typename std::enable_if<(sizeof...(Indices) == ndim) && (sizeof...(Indices) > 1),
                                        typename detail::return_type<true, T, Indices...>::type>::type
                operator()(Indices... idx) const
                {
                    return const_cast<marray_base&>(*this)(idx...);
                }

                pointer data()
                {
                    return mem_.data_;
                }

                const_pointer data() const
                {
                    return mem_.data_;
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                idx_type length() const
                {
                    return len_[0];
                }

                idx_type length(unsigned dim) const
                {
                    return len_[dim];
                }

                const std::array<idx_type, ndim>& lengths() const
                {
                    return len_;
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                stride_type stride() const
                {
                    return stride_[0];
                }

                stride_type stride(unsigned dim) const
                {
                    return stride_[dim];
                }

                const std::array<stride_type, ndim>& strides() const
                {
                    return stride_;
                }

                /* const view only */
                operator const marray_view<T,ndim>&() const { return static_cast<const marray_view<T,ndim>&>(*this); }

                void swap(marray&& other)
                {
                    swap(other);
                }

                void swap(marray& other)
                {
                    using std::swap;
                    swap(this->mem_.data_, other.mem_.data_);
                    swap(this->mem_.size_, other.mem_.size_);
                    swap(this->len_,       other.len_);
                    swap(this->stride_,    other.stride_);
                    swap(this->layout_,    other.layout_);
                }

                friend void swap(marray& a, marray& b)
                {
                    a.swap(b);
                }
        };
    }

    template <typename T, unsigned ndim>
    using const_marray_view = detail::marray_base<true, true, T, ndim>;

    template <typename T, unsigned ndim>
    using marray_view = detail::marray_base<false, true, T, ndim>;

    template <typename T, unsigned ndim, typename Allocator=aligned_allocator<T,MARRAY_BASE_ALIGNMENT>>
    using marray = detail::marray_base<false, false, T, ndim, Allocator>;

    template <typename T, int N> void copy(const marray_view<T,N>& a, marray_view<T,N>& b);

    /*
     * Convenient names for 1- and 2-dimensional array types.
     */
    template <typename T> using const_row_view = const_marray_view<T, 1>;
    template <typename T> using row_view = marray_view<T, 1>;
    template <typename T> using row = marray<T, 1>;

    template <typename T> using const_matrix_view = const_marray_view<T, 2>;
    template <typename T> using matrix_view = marray_view<T, 2>;
    template <typename T> using matrix = marray<T, 2>;

    template <typename T, int N>
    void copy(const marray_view<T,N>& a, marray_view<T,N>& b)
    {
        assert(a.lengths() == b.lengths());

        auto it = make_iterator(a.lengths(), a.strides(), b.strides());
        auto a_ = a.data();
        auto b_ = b.data();
        while (it.nextIteration(a_, b_)) *b_ = *a_;
    }

    template <typename T, int N>
    void copy(const marray_view<T,N>& a, marray_view<T,N>&& b)
    {
        copy(a, b);
    }

    template <typename T>
    matrix_view<T> operator^(matrix_view<T>& m, transpose_t t)
    {
        return matrix_view<T>(m.length(1), m.length(0), m.data(),
                              m.stride(1), m.stride(0));
    }

    template <typename T>
    matrix_view<T> operator^(matrix_view<T>&& m, transpose_t t)
    {
        return m^transpose::T;
    }

    template <typename T>
    const_matrix_view<T> operator^(const matrix_view<T>& m, transpose_t t)
    {
        return const_cast<matrix_view<T>&>(m)^transpose::T;
    }

    template <typename T>
    void gemm(const T& alpha, const matrix_view<T>& a,
                              const matrix_view<T>& b,
              const T&  beta,       matrix_view<T>& c)
    {
        //TODO
    }

    template <typename T>
    void gemm(const T& alpha, const matrix_view<T>&  a,
                              const matrix_view<T>&  b,
              const T&  beta,       matrix_view<T>&& c)
    {
        gemm(alpha, a, b, beta, c);
    }

    template <typename U>
    void gemm(char transa, char transb,
              const U& alpha, const matrix_view<U>& a,
                              const matrix_view<U>& b,
              const U&  beta,       matrix_view<U>& c)
    {
        using transpose::T;

        if (toupper(transa) == 'T')
        {
            if (toupper(transb) == 'T') gemm(alpha, a^T, b^T, beta, c);
            else                        gemm(alpha, a^T, b  , beta, c);
        }
        else
        {
            if (toupper(transb) == 'T') gemm(alpha, a  , b^T, beta, c);
            else                        gemm(alpha, a  , b  , beta, c);
        }
    }

    template <typename T>
    void gemm(char transa, char transb,
              const T& alpha, const matrix_view<T>&  a,
                              const matrix_view<T>&  b,
              const T&  beta,       matrix_view<T>&& c)
    {
        gemm(transa, transb, alpha, a, b, beta, c);
    }
}

#endif
