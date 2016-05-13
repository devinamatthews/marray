#ifndef _MARRAY_MARRAY_HPP_
#define _MARRAY_MARRAY_HPP_

#define MARRAY_TEST(...) break me!

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

    template <typename T, unsigned ndim> class const_marray_view;
    template <typename T, unsigned ndim> class marray_view;
    template <typename T, unsigned ndim> class marray;
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
        template <typename T, unsigned ndim> class marray_base;
        template <typename T_, unsigned ndim_, unsigned dim_> class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> class marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> class const_marray_slice;

        /*
         * Represents a part of an array, where the first dim-1 out of ndim
         * dimensions have been indexed into. This type may be implicity converted
         * to an array view or further indexed. This particular reference type
         * is explicitly const, and may only be used to read the original array.
         */
        template <typename T, unsigned ndim, unsigned dim>
        class const_marray_ref
        {
            template <typename T_, unsigned ndim_> friend class marray_base;
            template <typename T_, unsigned ndim_> friend class const_marray_view;
            template <typename T_, unsigned ndim_> friend class marray_view;
            template <typename T_, unsigned ndim_> friend class marray;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;

            protected:
                typedef marray_base<T, ndim> base;

                typedef typename base::stride_type stride_type;
                typedef typename base::idx_type idx_type;
                typedef typename base::value_type value_type;
                typedef typename base::pointer pointer;
                typedef typename base::const_pointer const_pointer;
                typedef typename base::reference reference;
                typedef typename base::const_reference const_reference;

                marray_base<T, ndim>& array;
                stride_type idx;

                const_marray_ref(const const_marray_ref& other) = default;

                const_marray_ref(const marray_base<T, ndim>& array, stride_type idx, idx_type i)
                : array(const_cast<marray_base<T, ndim>&>(array)), idx(idx+i*array.stride_[dim-2]) {}

                template <unsigned ndim_>
                marray_ref<T, ndim, dim>& operator=(const const_marray_ref<T, ndim_, ndim_-ndim+dim>& other)
                {
                    copy(other, *this);
                    return static_cast<marray_ref<T, ndim, dim>&>(*this);
                }

                template <unsigned ndim_, unsigned newdim_>
                marray_ref<T, ndim, dim>& operator=(const const_marray_slice<T, ndim_, ndim_-ndim+dim+newdim_, newdim_>& other)
                {
                    copy(other, *this);
                    return static_cast<marray_ref<T, ndim, dim>&>(*this);
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                marray_ref<T, ndim, dim+1> operator[](idx_type i)
                {
                    assert(i >= 0 && i < array.len_[dim-1]);
                    return {array, idx, i};
                }

                template <typename I, int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                marray_view<T, 1> operator[](const range_t<I>& x)
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[ndim-1]);
                    const_marray_view<T, 1> ret(*this);
                    ret.len_[0] = x.size();
                    ret.mem_.data_ += ret.stride_[0]*x.front();
                    return ret;
                }

                template <typename I, int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                marray_slice<T, ndim, dim+1, 1> operator[](const range_t<I>& x)
                {
                    return {array, idx, {}, {}, x};
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                marray_view<T, 1> operator[](const slice::all_t& x)
                {
                    return *this;
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                marray_slice<T, ndim, dim+1, 1> operator[](const slice::all_t& x) const
                {
                    return {array, idx, {}, {}, range(idx_type(), array.len_[dim-1])};
                }

                pointer data()
                {
                    return array.data_+idx;
                }

                operator marray_view<T, ndim-dim+1>()
                {
                    std::array<idx_type, ndim-dim+1> len;
                    std::array<stride_type, ndim-dim+1> stride;
                    std::copy_n(array.len_.begin(), ndim-dim+1, len.begin());
                    std::copy_n(array.stride_.begin(), ndim-dim+1, stride.begin());
                    return {len, data(), stride};
                }

            public:
                const_marray_ref& operator=(const const_marray_ref&) = delete;

                template <int diff=ndim-dim, typename=typename std::enable_if<diff==0>::type>
                const_reference operator[](idx_type i) const
                {
                    return const_cast<const_marray_ref&>(*this)[i];
                }

                template <int diff=ndim-dim, typename=typename std::enable_if<diff!=0>::type>
                const_marray_ref<T, ndim, dim+1> operator[](idx_type i) const
                {
                    return const_cast<const_marray_ref&>(*this)[i];
                }

                template <typename I, int diff=ndim-dim, typename=typename std::enable_if<diff==0>::type>
                const_marray_view<T, 1> operator[](const range_t<I>& x) const
                {
                    return const_cast<const_marray_ref&>(*this)[x];
                }

                template <typename I, int diff=ndim-dim, typename=typename std::enable_if<diff!=0>::type>
                const_marray_slice<T, ndim, dim+1, 1> operator[](const range_t<I>& x) const
                {
                    return const_cast<const_marray_ref&>(*this)[x];
                }

                template <int diff=ndim-dim, typename=typename std::enable_if<diff==0>::type>
                const_marray_view<T, 1> operator[](const slice::all_t& x) const
                {
                    return const_cast<const_marray_ref&>(*this)[x];
                }

                template <int diff=ndim-dim, typename=typename std::enable_if<diff!=0>::type>
                const_marray_slice<T, ndim, dim+1, 1> operator[](const slice::all_t& x) const
                {
                    return const_cast<const_marray_ref&>(*this)[x];
                }

                const_pointer data() const
                {
                    return const_cast<const_marray_ref&>(*this).data();
                }

                operator const_marray_view<T, ndim-dim+1>() const
                {
                    return operator marray_view<T, ndim-dim+1>();
                }
        };

        template <typename T, unsigned ndim, unsigned dim>
        class marray_ref : public const_marray_ref<T, ndim, dim>
        {
            template <typename T_, unsigned ndim_> friend class marray_base;
            template <typename T_, unsigned ndim_> friend class const_marray_view;
            template <typename T_, unsigned ndim_> friend class marray_view;
            template <typename T_, unsigned ndim_> friend class marray;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;

            protected:
                typedef marray_base<T, ndim> base;
                typedef const_marray_ref<T, ndim, dim> parent;

                typedef typename base::stride_type stride_type;
                typedef typename base::idx_type idx_type;
                typedef typename base::value_type value_type;
                typedef typename base::pointer pointer;
                typedef typename base::const_pointer const_pointer;
                typedef typename base::reference reference;
                typedef typename base::const_reference const_reference;

                using parent::array;
                using parent::idx;

                marray_ref(const marray_ref& other) = default;

                marray_ref(marray_base<T, ndim>& array, stride_type idx, idx_type i)
                : parent(array, idx, i) {}

            public:
                marray_ref& operator=(const marray_ref&) = delete;

                using parent::operator=;

                using parent::operator[];

                using parent::data;
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
            template <typename T_, unsigned ndim_> friend class marray_base;
            template <typename T_, unsigned ndim_> friend class const_marray_view;
            template <typename T_, unsigned ndim_> friend class marray_view;
            template <typename T_, unsigned ndim_> friend class marray;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;

            protected:
                typedef marray_base<T, ndim> base;

                typedef typename base::stride_type stride_type;
                typedef typename base::idx_type idx_type;
                typedef typename base::value_type value_type;
                typedef typename base::pointer pointer;
                typedef typename base::const_pointer const_pointer;
                typedef typename base::reference reference;
                typedef typename base::const_reference const_reference;

                marray_base<T, ndim>& array;
                stride_type idx;
                std::array<unsigned, newdim> dims;
                std::array<idx_type, newdim> lens;

                const_marray_slice(const const_marray_slice& other) = default;

                const_marray_slice(const marray_base<T, ndim>& array, stride_type idx,
                                   const std::array<unsigned,newdim>& dims,
                                   const std::array<idx_type,newdim>& lens, idx_type i)
                : array(marray_base<T, ndim>&(array)), idx(idx+i*array.stride[dim-2]), dims(dims), lens(lens) {}

                template <typename I>
                const_marray_slice(const marray_base<T, ndim>& array, stride_type idx,
                                   const std::array<unsigned,newdim-1>& dims,
                                   const std::array<idx_type,newdim-1>& lens, const range_t<I>& range_)
                : array(marray_base<T, ndim>&(array)), idx(idx+array.stride_[dim-2]*range_.front())
                {
                    std::copy(dims.begin(), dims.end(), this->dims.begin());
                    this->dims.back() = dim-2;
                    std::copy(lens.begin(), lens.end(), this->lens.begin());
                    this->lens.back() = range_.size();
                }

                template <unsigned ndim_>
                marray_slice<T, ndim, dim, newdim>& operator=(const const_marray_ref<T, ndim_, ndim_-ndim+dim-newdim>& other)
                {
                    copy(other, *this);
                    return static_cast<marray_slice<T, ndim, dim, newdim>&>(*this);
                }

                template <unsigned ndim_, unsigned newdim_>
                marray_slice& operator=(const const_marray_slice<T, ndim_, ndim_-ndim+dim-newdim+newdim_, newdim_>& other)
                {
                    copy(other, *this);
                    return static_cast<marray_slice<T, ndim, dim, newdim>&>(*this);
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                marray_view<T, newdim> operator[](idx_type i)
                {
                    assert(i >= 0 && i < array.len_[ndim-1]);
                    marray<T, newdim> ret;
                    ret.data_ = array.data_+idx+i*array.stride_[ndim-1];
                    ret.size_ = 1;
                    ret.is_view_ = true;
                    ret.layout_ = array.layout_;
                    for (unsigned dim = 0;dim < newdim;dim++)
                    {
                        ret.len_[dim] = lens[dim];
                        ret.stride_[dim] = array.stride_[dims[dim]];
                        ret.size_ *= ret.len_[dim];
                    }
                    return ret;
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                marray_slice<T, ndim, dim+1, newdim> operator[](idx_type i)
                {
                    assert(i >= 0 && i < array.len[dim-1]);
                    return {array, idx, dims, lens, i};
                }

                template <typename I, int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                marray_view<T, newdim+1> operator[](const range_t<I>& x)
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[ndim-1]);
                    marray<T, newdim+1> ret;
                    ret.data_ = array.data_+idx+array.stride_[ndim-1]*x.front();
                    ret.size_ = 1;
                    ret.is_view_ = true;
                    ret.layout_ = array.layout_;
                    for (unsigned dim = 0;dim < newdim;dim++)
                    {
                        ret.len_[dim] = lens[dim];
                        ret.stride_[dim] = array.stride_[dims[dim]];
                        ret.size_ *= ret.len_[dim];
                    }
                    ret.len_[newdim] = x.size();
                    ret.stride_[newdim] = array.stride_[ndim-1];
                    ret.size_ *= ret.len_[newdim];
                    return ret;
                }

                template <typename I, int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                marray_slice<T, ndim, dim+1, newdim+1> operator[](const range_t<I>& x)
                {
                    assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array.len_[dim-1]);
                    return {array, idx, dims, lens, x};
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                marray_view<T, newdim+1> operator[](const slice::all_t& x)
                {
                    return *this;
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                marray_slice<T, ndim, dim+1, newdim+1> operator[](const slice::all_t& x)
                {
                    return {array, idx, dims, lens, range(idx_type(), array.len_[dim-1])};
                }

                pointer data()
                {
                    return array.data_+idx;
                }

                operator marray_view<T, ndim+newdim-dim+1>()
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

            public:
                const_marray_slice& operator=(const const_marray_slice&) = delete;

                template <int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                const_marray_view<T, newdim> operator[](idx_type i) const
                {
                    return const_cast<const_marray_slice&>(*this)[i];
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                const_marray_slice<T, ndim, dim+1, newdim> operator[](idx_type i) const
                {
                    return const_cast<const_marray_slice&>(*this)[i];
                }

                template <typename I, int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                const_marray_view<T, newdim+1> operator[](const range_t<I>& x) const
                {
                    return const_cast<const_marray_slice&>(*this)[x];
                }

                template <typename I, int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                const_marray_slice<T, ndim, dim+1, newdim+1> operator[](const range_t<I>& x) const
                {
                    return const_cast<const_marray_slice&>(*this)[x];
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff==0>::type>
                const_marray_view<T, newdim+1> operator[](const slice::all_t& x) const
                {
                    return const_cast<const_marray_slice&>(*this)[x];
                }

                template <int diff=ndim-dim, typename=std::enable_if<diff!=0>::type>
                const_marray_slice<T, ndim, dim+1, newdim+1> operator[](const slice::all_t& x) const
                {
                    return const_cast<const_marray_slice&>(*this)[x];
                }

                const_pointer data() const
                {
                    return const_cast<const_marray_slice&>(*this).data();
                }

                operator const_marray_view<T, ndim+newdim-dim+1>() const
                {
                    return operator marray_view<T, ndim+newdim-dim+1>();
                }
        };

        template <typename T, unsigned ndim, unsigned dim, unsigned newdim>
        class marray_slice : public const_marray_slice<T, ndim, dim, newdim>
        {
            template <typename T_, unsigned ndim_> friend class marray_base;
            template <typename T_, unsigned ndim_> friend class const_marray_view;
            template <typename T_, unsigned ndim_> friend class marray_view;
            template <typename T_, unsigned ndim_> friend class marray;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;

            protected:
                typedef marray_base<T, ndim> base;
                typedef const_marray_slice<T, ndim, dim, newdim> parent;

                typedef typename base::stride_type stride_type;
                typedef typename base::idx_type idx_type;
                typedef typename base::value_type value_type;
                typedef typename base::pointer pointer;
                typedef typename base::const_pointer const_pointer;
                typedef typename base::reference reference;
                typedef typename base::const_reference const_reference;

                using parent::array;
                using parent::idx;
                using parent::dims;
                using parent::lens;

                marray_slice(const marray_slice& other) = default;

                marray_slice(marray_base<T, ndim>& array, stride_type idx,
                             const std::array<unsigned,newdim>& dims,
                             const std::array<idx_type,newdim>& lens, idx_type i)
                : parent(array, idx, dims, lens, i) {}

                template <typename I>
                marray_slice(marray_base<T, ndim>& array, stride_type idx,
                             const std::array<unsigned,newdim-1>& dims,
                             const std::array<idx_type,newdim-1>& lens, const range_t<I>& x)
                : parent(array, idx, dims, lens, x) {}

            public:
                marray_slice& operator=(const marray_slice&) = delete;

                using parent::operator=;

                using parent::operator[];

                using parent::data;
        };

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

        template <bool isConst, typename T, unsigned N>
        struct marray_type<false, T, N> { typedef marray_view<T, N> type; };
        
        template <typename T, unsigned N>
        struct marray_type<true, T, N> { typedef const_marray_view<T, N> type; };

        template <typename T>
        struct marray_type<false, T, 0> { typedef typename marray_base<T, 1>::reference type; };

        template <typename T>
        struct marray_type<true, T, 0> { typedef typename marray_base<T, 1>::const_reference type; };

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
        struct num_slices<Index, Indices>
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
        using are_integral = are_integral_helper<N, 1, Condition, Args...>;

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

        template <typename T, unsigned N, typename Tail, typename... Args>
        using set_len = set_len_helper<T, N, 1, Tail, Args...>;

        /*
         * Helper class to determine if the arguments Args... are of one of the
         * forms which are acceptable for passing to marray<T, ndim>::reset.
         */
        template <typename T, unsigned ndim>
        struct are_reset_args
        {
            template <typename... Args>
            struct direct_helper : std::false_type {};

            template <>
            struct direct_helper<> : std::true_type {};

            template <typename Arg>
            struct direct_helper<Arg>
            : std::integral_constant<bool, std::is_convertible<Arg,T>::value ||
                                           std::is_same<Arg,uninitialized_t>::value ||
                                           std::is_same<Arg,Layout>::value> {};

            template <typename Arg1, typename Arg2>
            struct direct_helper<Arg1, Arg2>
            : std::integral_constant<bool, (std::is_convertible<Arg1,T>::value ||
                                            std::is_same<Arg1,uninitialized_t>::value) &&
                                           std::is_same<Arg2,Layout>::value> {};

            template <typename... Args>
            using direct = std::integral_constant<bool,
                are_integral<ndim, direct_helper,
                             typename std::decay<Args>::type...>::value>;

            template <typename... Args>
            struct const_view_helper : std::false_type {};

            template <typename Arg>
            struct const_view_helper<Arg> : std::is_convertible<Arg,const T*> {};

            template <typename Arg>
            struct const_view_helper<Arg, Layout> : std::is_convertible<Arg,const T*> {};

            template <typename Arg, typename... Args>
            struct const_view_helper<Arg, Args...>
            : std::integral_constant<bool, std::is_convertible<Arg,const T*>::value &&
                                           are_integral<ndim, are_empty, Args...>::value> {};

            template <typename... Args>
            using const_view = std::integral_constant<bool,
                are_integral<ndim, const_view_helper,
                             typename std::decay<Args>::type...>::value>;

            template <typename... Args>
            struct view_helper : std::false_type {};

            template <typename Arg>
            struct view_helper<Arg> : std::is_convertible<Arg,T*> {};

            template <typename Arg>
            struct view_helper<Arg, Layout> : std::is_convertible<Arg,T*> {};

            template <typename Arg, typename... Args>
            struct view_helper<Arg, Args...>
            : std::integral_constant<bool, std::is_convertible<Arg,T*>::value &&
                                           are_integral<ndim, are_empty, Args...>::value> {};

            template <typename... Args>
            using view = std::integral_constant<bool,
                are_integral<ndim, view_helper,
                             typename std::decay<Args>::type...>::value>;
        };

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

        template <typename T, typename ndim, typename... Args>
        void reset(const_marray_view<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<T, ndim, reset_helper, Args...> tmp(array, len, args...);
        }

        template <typename T, typename ndim, typename... Args>
        void reset(marray_view<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<T, ndim, reset_helper, Args...> tmp(array, len, args...);
        }

        template <typename T, typename ndim, typename... Args>
        void reset(marray<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<T, ndim, reset_helper, Args...> tmp(array, len, args...);
        }

        /*
         * Helper class to determine if the arguments Args... are of one of the
         * forms which are acceptable for passing to marray<T, ndim>::resize.
         */
        template <typename T, unsigned ndim>
        struct are_resize_args
        {
            template <typename... Args>
            struct direct_helper : std::false_type {};

            template <>
            struct direct_helper<> : std::true_type {};

            template <typename Arg>
            struct direct_helper<Arg>
            : std::integral_constant<bool, std::is_convertible<Arg,T>::value> {};

            template <typename... Args>
            using direct =  std::integral_constant<bool,
                 are_integral<ndim, direct_helper,
                              typename std::decay<Args>::type...>::value>;
        };

        struct resize_helper
        {
            template <typename T, unsigned ndim>
            resize_helper(marray<T, ndim>& array, std::array<size_t, ndim>& len, const T& val=T())
            {
                array.resize(len, val);
            }
        };

        template <typename T, typename ndim, typename... Args>
        void resize(marray<T, ndim>& array, const Args&... args)
        {
            std::array<size_t, ndim> len;
            set_len<T, ndim, resize_helper, Args...> tmp(array, len, args...);
        }

        template <typename T, unsigned ndim, typename Allocator=aligned_allocator<T,MARRAY_BASE_ALIGNMENT>>
        class marray_base
        {
            static_assert(ndim > 0, "0-dimensional marrays are not allowed.");

            template <typename T_, unsigned ndim_> friend class detail::marray_base;
            template <typename T_, unsigned ndim_> friend class const_marray_view;
            template <typename T_, unsigned ndim_> friend class marray_view;
            template <typename T_, unsigned ndim_> friend class marray;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
            template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
            template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
            template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

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

                    mem_t(Allocator alloc) : Allocator(std::move(alloc)) {}

                    ~mem_t() { clear(); }

                    void allocate(size_type size, unitialized_t u)
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
                        size_ = size;
                        is_alloced_ = true;
                        data_ = alloc.allocate(other.size_);
                    }

                    void clear()
                    {
                        if (is_alloced_ && data_)
                        {
                            for (size_type i = 0;i < size;i++)
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

                marray_base() {}

                marray_base(const marray_base& other) = delete;

                marray_base(marray_base&& other) = delete;

                marray_base& operator=(const marray_base& other) = delete;

                marray_base& operator=(marray_base&& other) = delete;

                marray_base& operator=(const marray_view<T,ndim>& other)
                {
                    copy(other, static_cast<marray_view<T,ndim>&>(*this));
                    return *this;
                }

                void view(const marray_base& other)
                {
                    mem_.set(other.mem_.data_);
                    len_ = other.len_;
                    stride_ = other.stride_;
                }

                void copy(const marray_base& other, Layout layout=DEFAULT)
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

                void copy(marray_base&& other)
                {
                    other.mem_.set(std::move(other.mem_));
                    len_ = other.len_;
                    stride_ = other.stride_;
                }

                template <typename U>
                void reset(const std::array<U, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
                {
                    reset(len, uninitialized, layout);
                    std::uninitialized_fill_n(mem_.data_, mem_.size_, val);
                }

                template <typename U>
                void reset(const std::array<U, ndim>& len, uninitialized_t u, Layout layout = DEFAULT)
                {
                    std::copy_n(len.begin(), ndim, len_.begin());

                    size_t size;
                    if (layout == ROW_MAJOR)
                    {
                        stride_[ndim-1] = 1;
                        for (unsigned i = ndim-1;i > 0;i--)
                        {
                            stride_[i-1] = stride_[i]*len_[i];
                        }
                        size = stride_[0]*len_[0];
                    }
                    else
                    {
                        stride_[0] = 1;
                        for (unsigned i = 1;i < ndim;i++)
                        {
                            stride_[i] = stride_[i-1]*len_[i-1];
                        }
                        size = stride_[ndim-1]*len_[ndim-1];
                    }

                    mem_.allocate(size, u);
                }

                template <typename U>
                marray_view<T, ndim> permute(const std::array<U, ndim>& perm)
                {
                    marray_view<T, ndim> view;

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
                typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, marray_view<T, ndim>>::type
                permute(Args&&... args)
                {
                    return permute(make_array(std::forward<Args>(args)...));
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                reference front()
                {
                    assert(length() > 0);
                    return mem_.data_[0];
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
                reference front(unsigned dim)
                {
                    assert(dim == 0);
                    return front();
                }

                template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
                marray_view<T, ndim-1> front(unsigned dim)
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
                reference operator()(idx_type i)
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

                pointer data()
                {
                    return mem_.data_;
                }

                void rotate_dim(unsigned dim, idx_type shift)
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
                    Iterator<idx_type, stride_type> it(sublen, substride);
                    while (it.nextIteration(p))
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

                template <typename U>
                void rotate(const std::array<U, ndim>& shift)
                {
                    for (unsigned dim = 0;dim < ndim;dim++)
                    {
                        rotate_dim(dim, shift[dim]);
                    }
                }

                template <typename... Args>
                typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value>::type
                rotate(Args&&... args)
                {
                    rotate(make_array(std::forward<Args>(args)...));
                }

            public:
                void reset()
                {
                    mem_.clear();
                    len_.fill(0);
                    stride_.fill(0);
                }

                template <typename U>
                void reset(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
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
                }

                template <typename U, typename V>
                void reset(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
                {
                    mem_.set(ptr);
                    std::copy_n(len.begin(), ndim, len_.begin());
                    std::copy_n(stride.begin(), ndim, stride_.begin());
                }

                template <typename U>
                const_marray_view<T, ndim> permute(const std::array<U, ndim>& perm) const
                {
                    return const_cast<marray_base&>(*this).permute(perm);
                }

                template <typename... Args>
                typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, const_marray_view<T, ndim>>::type
                permute(Args&&... args) const
                {
                    return permute(make_array(std::forward<Args>(args)...));
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
                const_reference operator()(idx_type i) const
                {
                    return (*this)[i];
                }

                template <typename... Indices>
                typename std::enable_if<(sizeof...(Indices) == ndim) && (sizeof...(Indices) > 1),
                                        typename detail::return_type<true, T, Indices...>::type>::type
                operator()(Indices... idx) const
                {
                    return const_cast<marray_base&>(*this)(idx...);
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
        };
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
    class const_marray_view : detail::marray_base<T,ndim>
    {
        static_assert(ndim > 0, "0-dimensional marrays are not allowed.");

        template <typename T_, unsigned ndim_> friend class detail::marray_base;
        template <typename T_, unsigned ndim_> friend class const_marray_view;
        template <typename T_, unsigned ndim_> friend class marray_view;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

        protected:
            typedef detail::marray_base<T,ndim> base;

        public:
            using typename base::idx_type;
            using typename base::stride_type;
            using typename base::value_type;
            using typename base::pointer;
            using typename base::const_pointer;
            using typename base::reference;
            using typename base::const_reference;

            const_marray_view() {}

            const_marray_view(const const_marray_view& other)
            {
                base::reset(other);
            }

            template <typename U>
            const_marray_view(const std::array<U, ndim>& len, const_pointer ptr, Layout layout = DEFAULT)
            {
                reset(len, ptr, layout);
            }

            template <typename U, typename V>
            const_marray_view(const std::array<U, ndim>& len, const_pointer ptr, const std::array<V, ndim>& stride)
            {
                reset(len, ptr, stride);
            }

            template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim, Args...>::value>::type>
            explicit const_marray_view(Args&&... args)
            {
                reset(std::forward<Args>(args)...);
            }

            template <typename U>
            void reset(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
            {
                base::reset(len, const_cast<pointer>(ptr), layout);
            }

            template <typename U, typename V>
            void reset(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
            {
                base::reset(len, const_cast<pointer>(ptr), stride);
            }

            template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim>::const_view<Args...>::value>::type>
            void reset(Args&&... args)
            {
                detail::reset(*this, std::forward<Args>(args)...);
            }

            operator const marray_view<T,ndim>&() const { return static_cast<const marray_view<T,ndim>&>(*this); }
    };

    template <typename T, unsigned ndim>
    class marray_view : public const_marray_view<T, ndim>
    {
        template <typename T_, unsigned ndim_> friend class detail::marray_base;
        template <typename T_, unsigned ndim_> friend class const_marray_view;
        template <typename T_, unsigned ndim_> friend class marray_view;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

        protected:
            typedef detail::marray_base<T,ndim> base;

        public:
            using typename base::idx_type;
            using typename base::stride_type;
            using typename base::value_type;
            using typename base::pointer;
            using typename base::const_pointer;
            using typename base::reference;
            using typename base::const_reference;

        public:
            marray_view() {}

            marray_view(marray_view& other)
            : const_marray_view<T, ndim>(other) {}

            marray_view(marray_view&& other)
            : const_marray_view<T, ndim>(other) {}

            template <typename U>
            marray_view(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
            : const_marray_view<T, ndim>(len, ptr, layout) {}

            template <typename U, typename V>
            marray_view(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
            : const_marray_view<T, ndim>(len, ptr, stride) {}

            template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim, Args...>::value>::type>
            marray_view(Args&&... args)
            : const_marray_view<T, ndim>(std::forward<Args>(args)...) {}

            template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim>::view<Args...>::value>::type>
            void reset(Args&&... args)
            {
                detail::reset(*this, std::forward<Args>(args)...);
            }

            using base::operator=;

            using base::permute;

            using base::front;

            using base::back;

            using base::operator[];

            using base::operator();

            using base::data;

            using base::rotate_dim;

            using base::rotate;
    };

    template <typename T, unsigned ndim>
    class marray : public marray_view<T, ndim>
    {
        template <typename T_, unsigned ndim_> friend class detail::marray_base;
        template <typename T_, unsigned ndim_> friend class const_marray_view;
        template <typename T_, unsigned ndim_> friend class marray_view;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
        template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

        protected:
            typedef detail::marray_base<T,ndim> base;

            Layout layout_;

        public:
            using typename base::idx_type;
            using typename base::stride_type;
            using typename base::value_type;
            using typename base::pointer;
            using typename base::const_pointer;
            using typename base::reference;
            using typename base::const_reference;

        public:
            marray() : layout_(DEFAULT) {}

            marray(const marray_view<T, ndim>& other, Layout layout=DEFAULT)
            : layout_(layout)
            {
                base::copy(other, layout);
            }

            marray(const marray& other, Layout layout=DEFAULT)
            : layout_(layout)
            {
                base::copy(other, layout);
            }

            marray(marray&& other, Layout layout=DEFAULT)
            : layout_(layout)
            {
                base::copy(std::move(other), layout);
            }

            template <typename U>
            marray(const std::array<U, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
            : layout_(layout)
            {
                base::reset(len, val, layout);
            }

            template <typename U>
            marray(const std::array<U, ndim>& len, uninitialized_t u, Layout layout = DEFAULT)
            : layout_(layout)
            {
                base::reset(len, u, layout);
            }

            template <typename... Args, typename std::enable_if<detail::are_reset_args<T, ndim, Args...>::value>::type>
            marray(Args&&... args)
            {
                reset(std::forward<Args>(args)...);
                layout_ = stride_[0] == 1 ? COLUMN_MAJOR : ROW_MAJOR;
            }

            using base::reset;

            void clear()
            {
                reset();
            }

            template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim>::direct<Args...>::value>::type>
            void reset(Args&&... args)
            {
                detail::reset(*this, std::forward<Args>(args)...);
            }

            template <typename U>
            void resize(const std::array<U, ndim>& len, const T& val=T())
            {
                marray a(std::move(*this));
                reset(len, val, layout_);
                marray_view<T,ndim> b(*this);

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

            template <typename... Args, typename=typename std::enable_if<detail::are_resize_args<T, ndim>::direct<Args...>::value>::type>
            void resize(Args&&... args)
            {
                detail::resize(*this, std::forward<Args>(args)...);
            }

            template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
            void push_back(const T& x)
            {
                resize(length()+1);
                back() = x;
            }

            template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
            void pop_back()
            {
                resize(length()-1);
            }

            void push_back(unsigned dim, const marray_view<T, ndim-1>& x)
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
