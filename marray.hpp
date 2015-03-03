#ifndef _MARRAY_MARRAY_HPP_
#define _MARRAY_MARRAY_HPP_

#include <vector>
#include <type_traits>
#include <array>
#include <cassert>
#include <algorithm>

#include "iterator.hpp"

namespace MArray
{

    namespace detail
    {
        template <typename T, typename... Args>
        struct first_type
        {
            typedef typename std::decay<T>::type type;
        };

        template <typename T>
        void consume_args(std::vector<T>& v) {}

        template <typename T, typename... Args>
        void consume_args(std::vector<T>& v, const T& t, const Args&... args)
        {
            v.push_back(t);
            consume_args(v, args...);
        }
    }

    /*
     * Create a vector from the specified elements, where the type of the vector
     * is taken from the first element.
     */
    template <typename... T> std::vector<typename detail::first_type<T...>::type>
    vec(const T&... a)
    {
        std::vector<typename detail::first_type<T...>::type> v;
        v.reserve(sizeof...(T));
        detail::consume_args(v, a...);
        return v;
    }

    template <typename T, unsigned ndim> class const_marray;
    template <typename T, unsigned ndim> class marray;
    template <typename T> class matrix;
    template <typename T> class row;
    template <typename T, unsigned ndim, unsigned dim> class const_marray_ref;
    template <typename T, unsigned ndim, unsigned dim> class marray_ref;
    template <typename T, unsigned ndim, unsigned dim, unsigned newdim> class const_marray_slice;
    template <typename T, unsigned ndim, unsigned dim, unsigned newdim> class marray_slice;

    namespace slice
    {
        struct slice_all_t {};
        constexpr slice_all_t all = slice_all_t();
    }

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

            const_marray_slice<T, ndim, dim+1, 1> operator[](const slice::slice_all_t& x) const
            {
                return const_marray_slice<T, ndim, dim+1, 1>(array_, idx, {}, dim);
            }
    };

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

            marray_slice<T, ndim, dim+1, 1> operator[](const slice::slice_all_t& x)
            {
                return marray_slice<T, ndim, dim+1, 1>(array_, idx, {}, dim);
            }
    };

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

            const_marray<T, 1> operator[](const slice::slice_all_t& x) const
            {
                return const_marray<T, 1>(*this);
            }
    };

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

            marray<T, 1> operator[](const slice::slice_all_t& x) const
            {
                return marray<T, 1>(*this);
            }
    };

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
            std::array<idx_type, newdim> dims;

            const_marray_slice(const const_marray_slice& other) = default;

            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               idx_type i, const std::array<idx_type,newdim-1>& dims)
            : array_(array_), idx(idx+i*array_.stride[dim-2]), dims(dims) {}

            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               const std::array<idx_type,newdim-1>& dims, unsigned dim_)
            : array_(array_), idx(idx)
            {
                std::copy(dims.begin(), dims.end(), this->dims.begin());
                this->dims.back() = dim_-1;
            }

        public:
            const_marray_slice& operator=(const const_marray_slice& other) = delete;

            const_marray_slice<T, ndim, dim+1, newdim> operator[](idx_type i) const
            {
                assert(i >= 0 && i < array_.len[dim-1]);
                return const_marray_slice<T, ndim, dim+1, newdim>(array_, idx, i, dims);
            }

            const_marray_slice<T, ndim, dim+1, newdim+1> operator[](const slice::slice_all_t& x) const
            {
                return const_marray_slice<T, ndim, dim+1, newdim+1>(array_, idx, dims, dim);
            }
    };

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

            marray_slice(const marray_slice& other) = default;

            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               idx_type i, const std::array<idx_type,newdim-1>& dims)
            : const_marray_slice<T, ndim, dim, newdim>(array_, idx, i, dims) {}

            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               const std::array<idx_type,newdim-1>& dims, unsigned dim_)
            : const_marray_slice<T, ndim, dim, newdim>(array_, idx, dims, dim_) {}

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
                assert(i >= 0 && i < array_.len[dim-1]);
                return marray_slice<T, ndim, dim+1, newdim>(array_, idx, i, dims);
            }

            marray_slice<T, ndim, dim+1, newdim+1> operator[](const slice::slice_all_t& x)
            {
                return marray_slice<T, ndim, dim+1, newdim+1>(array_, idx, dims, dim);
            }
    };

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
            std::array<idx_type, newdim> dims;

            const_marray_slice(const const_marray_slice& other) = default;

            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               idx_type i, const std::array<idx_type,newdim-1>& dims)
            : array_(array_), idx(idx+i*array_.stride[ndim-2]), dims(dims) {}

            const_marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                               const std::array<idx_type,newdim-1>& dims, unsigned dim)
            : array_(array_), idx(idx)
            {
                std::copy(dims.begin(), dims.end(), this->dims.begin());
                this->dims.back() = dim-1;
            }

        public:
            const_marray_slice& operator=(const const_marray_slice& other) = delete;

            const_marray<T, newdim> operator[](idx_type i) const
            {
                assert(i >= 0 && i < array_.len[ndim-1]);
                const_marray<T, newdim> ret;
                ret.data_ = array_.data_+idx+i*array_.stride[ndim-1];
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = array_.len_[dims[dim]];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                return ret;
            }

            const_marray<T, newdim+1> operator[](const slice::slice_all_t& x) const
            {
                const_marray<T, newdim+1> ret;
                ret.data_ = array_.data_+idx;
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = array_.len_[dims[dim]];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                ret.len_[newdim] = array_.len_[ndim-1];
                ret.stride_[newdim] = array_.stride_[ndim-1];
                ret.size_ *= ret.len_[newdim];
                return ret;
            }
    };

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

            marray_slice(const marray_slice& other) = default;

            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                         idx_type i, const std::array<idx_type,newdim-1>& dims)
            : const_marray_slice<T, ndim, ndim, newdim>(array_, idx, i, dims) {}

            marray_slice(const const_marray<T, ndim>& array_, size_type idx,
                         const std::array<idx_type,newdim-1>& dims, unsigned dim)
            : const_marray_slice<T, ndim, ndim, newdim>(array_, idx, dims, dim) {}

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
                assert(i >= 0 && i < array_.len[ndim-1]);
                marray<T, newdim> ret;
                ret.data_ = array_.data_+idx+i*array_.stride[ndim-1];
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = array_.len_[dims[dim]];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                return ret;
            }

            marray<T, newdim+1> operator[](const slice::slice_all_t& x)
            {
                marray<T, newdim+1> ret;
                ret.data_ = array_.data_+idx;
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = array_.len_[dims[dim]];
                    ret.stride_[dim] = array_.stride_[dims[dim]];
                    ret.size_ *= ret.len_[dim];
                }
                ret.len_[newdim] = array_.len_[ndim-1];
                ret.stride_[newdim] = array_.stride_[ndim-1];
                ret.size_ *= ret.len_[newdim];
                return ret;
            }
    };

    struct construct_view_t {};
    constexpr construct_view_t construct_view = construct_view_t();

    enum Layout {COLUMN_MAJOR, ROW_MAJOR, DEFAULT};
    constexpr Layout defaultLayout = ROW_MAJOR;

    namespace detail
    {
        template <bool isConst, typename T, unsigned N>
        struct marray_type { typedef marray<T, N> type; };
        
        template <typename T, unsigned N>
        struct marray_type<true, T, N> { typedef const_marray<T, N> type; };

        template <typename T>
        struct marray_type<false, T, 0> { typedef typename marray<T, 1>::reference type; };

        template <typename T>
        struct marray_type<true, T, 0> { typedef typename marray<T, 1>::const_reference type; };

        template <typename Index>
        struct is_slice { static constexpr unsigned N = 0; };

        template <>
        struct is_slice<slice::slice_all_t> { static constexpr unsigned N = 1; };

        template <typename T, typename Index, typename... Indices>
        struct const_return_type
        {
            static constexpr unsigned N = const_return_type<T, Indices...>::N + is_slice<Index>::N;
            typedef typename marray_type<true, T, N>::type type;
        };

        template <typename T, typename Index>
        struct const_return_type<T, Index>
        {
            static constexpr unsigned N = is_slice<Index>::N;
            typedef typename marray_type<true, T, N>::type type;
        };

        template <typename T, typename Index, typename... Indices>
        struct return_type
        {
            static constexpr unsigned N = return_type<T, Indices...>::N + is_slice<Index>::N;
            typedef typename marray_type<false, T, N>::type type;
        };

        template <typename T, typename Index>
        struct return_type<T, Index>
        {
            static constexpr unsigned N = is_slice<Index>::N;
            typedef typename marray_type<false, T, N>::type type;
        };

        template <typename RT, typename Index, typename... Indices>
        struct get_const_slice
        {
            template <typename ArrayOrSlice>
            RT operator()(const ArrayOrSlice& x, const Index& idx0, const Indices&... idx)
            {
                return get_const_slice<RT, Indices...>()(x[idx0], idx...);
            }
        };

        template <typename RT, typename Index>
        struct get_const_slice<RT, Index>
        {
            template <typename ArrayOrSlice>
            RT operator()(const ArrayOrSlice& x, const Index& idx0)
            {
                return x[idx0];
            }
        };

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
    }

    template <typename T, unsigned ndim>
    class const_marray
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

        public:
            typedef unsigned idx_type;
            typedef size_t size_type;
            typedef T value_type;
            typedef T* pointer;
            typedef const T* const_pointer;
            typedef T& reference;
            typedef const T& const_reference;

        protected:
            pointer data_;
            size_type size_;
            std::array<idx_type,ndim> len_;
            std::array<size_type,ndim> stride_;
            bool is_view_;
            Layout layout_;

            void reset(const std::array<idx_type, ndim>& len = std::array<idx_type, ndim>(),
                       const T& val = T(), Layout layout = DEFAULT)
            {
                size_type old_size = (is_view_ ? 0 : size_);
                len_ = len;

                if (layout == DEFAULT) layout = defaultLayout;
                if (layout == ROW_MAJOR)
                {
                    stride_[ndim-1] = 1;
                    for (unsigned i = ndim-1;i > 0;i--)
                    {
                        stride_[i-1] = stride_[i]*len[i];
                    }
                    size_ = stride_[0]*len[0];
                }
                else
                {
                    stride_[0] = 1;
                    for (unsigned i = 1;i < ndim;i++)
                    {
                        stride_[i] = stride_[i-1]*len[i-1];
                    }
                    size_ = stride_[ndim-1]*len[ndim-1];
                }

                if (size_ > old_size || is_view_)
                {
                    if (!is_view_ && data_) delete[] data_;
                    data_ = new T[size_];
                }
                std::fill(data_, data_+size_, val);

                is_view_ = false;
            }

        public:
            const_marray(const const_marray& other)
            : data_(NULL), size_(other.size_), len_(other.len_),
              stride_(other.stride_), is_view_(other.is_view_),
              layout_(other.layout_)
            {
                if (is_view_)
                {
                    data_ = other.data_;
                }
                else if (size_ > 0)
                {
                    data_ = new T[size_];
                    std::copy(other.data_, other.data_+size_, data_);
                }
            }

            const_marray(const_marray&& other)
            : data_(other.data_), size_(other.size_), len_(other.len_),
              stride_(other.stride_), is_view_(other.is_view_),
              layout_(other.layout_)
            {
                if (!is_view_)
                {
                    other.data_ = NULL;
                    other.reset();
                }
            }

            const_marray(const const_marray& other, const construct_view_t& cv)
            : data_(other.data_), size_(other.size_), len_(other.len_),
              stride_(other.stride_), is_view_(true),
              layout_(other.layout_) {}

            explicit const_marray(const std::array<idx_type, ndim>& len = std::array<idx_type, ndim>(),
                                  const T& val = T(), Layout layout = DEFAULT)
            : data_(NULL), len_(len), stride_(), is_view_(false), layout_(layout)
            {
                if (layout_ == DEFAULT) layout_ = defaultLayout;
                if (layout_ == ROW_MAJOR)
                {
                    stride_[ndim-1] = 1;
                    for (unsigned i = ndim-1;i > 0;i--)
                    {
                        stride_[i-1] = stride_[i]*len[i];
                    }
                    size_ = stride_[0]*len[0];
                }
                else
                {
                    stride_[0] = 1;
                    for (unsigned i = 1;i < ndim;i++)
                    {
                        stride_[i] = stride_[i-1]*len[i-1];
                    }
                    size_ = stride_[ndim-1]*len[ndim-1];
                }

                if (size_ > 0)
                {
                    data_ = new T[size_];
                    std::fill(data_, data_+size_, val);
                }
            }

            template <unsigned ndim_>
            const_marray(const const_marray_ref<T, ndim_, ndim_-ndim+1>& other)
            : data_(other.array_.data_+other.idx), layout_(other.array_.layout_),
              is_view_(true)
            {
                std::copy(other.array_.len_.begin()+ndim_-ndim, other.array_.len_.end(), len_.begin());
                std::copy(other.array_.stride_.begin()+ndim_-ndim, other.array_.stride_.end(), stride_.begin());

                size_ = 1;
                for (idx_type i = 0;i < ndim;i++) size_ *= len_[i];
            }

            template <unsigned ndim_, unsigned newdim_>
            const_marray(const const_marray_slice<T, ndim_, ndim_-ndim+1+newdim_, newdim_>& other)
            : data_(other.array_.data_+other.idx), layout_(other.array_.layout_),
              is_view_(true)
            {
                for (unsigned i = 0;i < newdim_;i++)
                {
                    len_[i] = other.array_.len_[other.dims[i]];
                    stride_[i] = other.array_.stride_[other.dims[i]];
                }
                std::copy(other.array_.len_.begin()+ndim_-ndim+newdim_, other.array_.len_.end(), len_.begin()+newdim_);
                std::copy(other.array_.stride_.begin()+ndim_-ndim+newdim_, other.array_.stride_.end(), stride_.begin()+newdim_);

                size_ = 1;
                for (idx_type i = 0;i < ndim;i++) size_ *= len_[i];
            }

            ~const_marray()
            {
                if (!is_view_ && data_) delete[] data_;
            }

            bool isView() const
            {
                return is_view_;
            }

            const_marray<T, ndim-1> front(unsigned dim) const
            {
                assert(dim >= 0 && dim < ndim);
                assert(len_[dim] > 0);

                const_marray<T, ndim-1> view;
                view.data_ = data_;
                view.stride_ = stride_;
                view.stride_.erase(view.stride_.begin()+dim);
                view.len_ = len_;
                view.len_.erase(view.len_.begin()+dim);
                view.is_view_ = true;
                view.layout_ = layout_;
                view.size_ = size_/len_[ndim];

                return view;
            }

            const_marray<T, ndim-1> back(unsigned dim) const
            {
                assert(dim >= 0 && dim < ndim);
                assert(len_[dim] > 0);

                const_marray<T, ndim-1> view;
                view.data_ = data_+(len_[dim]-1)*stride_[dim];
                view.stride_ = stride_;
                view.stride_.erase(view.stride_.begin()+dim);
                view.len_ = len_;
                view.len_.erase(view.len_.begin()+dim);
                view.is_view_ = true;
                view.layout_ = layout_;
                view.size_ = size_/len_[ndim];

                return view;
            }

            const_marray& operator=(const const_marray& other) = delete;

            const_marray_ref<T, ndim, 2> operator[](idx_type i) const
            {
                assert(i >= 0 && i < len_[0]);
                return const_marray_ref<T, ndim, 2>(*this, (size_type)0, i);
            }

            template <typename... Indices>
            typename std::enable_if<sizeof...(Indices) == ndim, typename detail::const_return_type<T, Indices...>::type>::type
            operator()(const Indices&... idx) const
            {
                return detail::get_const_slice<typename detail::const_return_type<T, Indices...>::type, Indices...>()(*this, idx...);
            }

            const_pointer data() const
            {
                return data_;
            }

            const idx_type& length(unsigned dim) const
            {
                return len_[dim];
            }

            const std::vector<idx_type>& lengths() const
            {
                return len_;
            }

            const size_type& stride(unsigned dim) const
            {
                return stride_[dim];
            }

            const std::vector<size_type>& strides() const
            {
                return stride_;
            }

            size_type size() const
            {
                return size_;
            }
    };

    template <typename T, unsigned ndim>
    class marray : public const_marray<T, ndim>
    {
        template <typename T_, unsigned ndim_> friend class const_marray;
        template <typename T_, unsigned ndim_> friend class marray;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
        template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;

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

            marray(const const_marray<T, ndim>& other)
            : const_marray<T, ndim>(other) {}

        public:
            marray(marray& other)
            : const_marray<T, ndim>(other) {}

            /*
             * Don't allow for now.
             */
            //marray(const const_marray<T, ndim>& other)
            //: const_marray<T, ndim>(other)
            //{
            //    unView();
            //}

            marray(marray&& other)
            : const_marray<T, ndim>(move(other)) {}

            marray(marray& other, const construct_view_t& cv)
            : const_marray<T, ndim>(other, cv) {}

            explicit marray(const std::array<idx_type, ndim>& len = std::array<idx_type, ndim>(),
                           const T& val = T(), Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, val, layout) {}

            template <unsigned ndim_>
            marray(const marray_ref<T, ndim_, ndim_-ndim+1>& other)
            : const_marray<T, ndim>(other) {}

            template <unsigned ndim_, unsigned newdim_>
            marray(const marray_slice<T, ndim_, ndim_-ndim+1+newdim_, newdim_>& other)
            : const_marray<T, ndim>(other) {}

            using const_marray<T, ndim>::reset;

            void clear()
            {
                reset();
            }

            void unView()
            {
                if (!is_view_) return;

                const_marray<T, ndim> old(*this);
                reset(len_, T(), layout_);
                copy(old, *this);
            }

            void push_back(unsigned dim, const const_marray<T, ndim-1>& x)
            {
                assert(dim >= 0 && dim < ndim);

                for (unsigned i = 0, j = 0;i < ndim;i++)
                {
                    if (i != dim)
                    {
                        assert(len_[i] == x.len_[j++]);
                    }
                }

                const_marray<T, ndim> old(*this);

                len_[dim]++;
                reset(len_, T(), layout_);

                T* old_data = old.data();
                T* new_data = data();

                for (Iterator<idx_type, size_type> it(old.len_, old.stride_, stride_);
                     it.nextIteration(old_data, new_data);) *new_data = *old_data;

                T* old_slice = x.data();
                T* new_slice = new_data;

                for (unsigned i = 0;i < ndim;i++)
                {
                    new_slice += (old.len_[i]-1)*stride_[i];
                }

                std::vector<size_type> slice_stride_ = stride_;
                slice_stride_.erase(slice_stride_.begin()+dim);

                for (Iterator<idx_type, size_type> it(x.len_, x.stride_, slice_stride_);
                     it.nextIteration(old_slice, new_slice);) *new_slice = *old_slice;
            }

            void pop_back(unsigned dim)
            {
                assert(dim >= 0 && dim < ndim);
                assert(len_[dim] > 0);

                const_marray<T, ndim> old(*this);

                len_[dim]--;
                reset(len_, T(), layout_);

                T* old_data = old.data();
                T* new_data = data();

                for (Iterator<idx_type, size_type> it(len_, old.stride_, stride_);
                     it.nextIteration(old_data, new_data);) *new_data = *old_data;
            }

            marray<T, ndim-1> front(unsigned dim)
            {
                return const_marray<T, ndim>::front(dim);
            }

            marray<T, ndim-1> back(unsigned dim)
            {
                return const_marray<T, ndim>::front(dim);
            }

            marray& operator=(const marray& other)
            {
                copy(other, *this);
                return *this;
            }

            marray& operator=(const const_marray<T, ndim>& other)
            {
                copy(other, *this);
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

            using const_marray<T, ndim>::operator[];

            marray_ref<T, ndim, 2> operator[](idx_type i)
            {
                assert(i >= 0 && i < len_[0]);
                return marray_ref<T, ndim, 2>(*this, (size_type)0, i);
            }

            using const_marray<T, ndim>::operator();

            template <typename... Indices>
            typename std::enable_if<sizeof...(Indices) == ndim, typename detail::return_type<T, Indices...>::type>::type
            operator()(const Indices&... idx)
            {
                return detail::get_slice<typename detail::return_type<T, Indices...>::type, Indices...>()(*this, idx...);
            }

            pointer data()
            {
                return data_;
            }

            friend void copy(const const_marray<T, ndim>& a, marray&& b)
            {
                copy(a, b);
            }

            friend void copy(const const_marray<T, ndim>& a, marray& b)
            {
                assert(a.len_ == b.len_);

                if (!a.is_view_ && !b.is_view_ && a.layout_ == b.layout_)
                {
                    std::copy(a.data_, a.data_+a.size_, b.data_);
                }
                else
                {
                    Iterator<idx_type, size_type> it(a.len_, a.stride_, b.stride_);
                    const_pointer a_ = a.data();
                          pointer b_ = b.data();
                    while (it.nextIteration(a_, b_)) *b_ = *a_;
                }
            }

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
            pointer data_;
            idx_type size_;
            size_type stride_;
            bool is_view_;

            void reset(idx_type n = 0, const T& val = T())
            {
                idx_type old_size = (is_view_ ? 0 : size_);
                size_ = n;

                if (size_ > old_size)
                {
                    if (!is_view_ && data_) delete[] data_;
                    data_ = new T[size_];
                }
                std::fill(data_, data_+size_, val);

                is_view_ = false;
            }

        public:
            const_marray(const const_marray& other)
            : data_(NULL), size_(other.size_), stride_(other.stride_),
              is_view_(other.is_view_)
            {
                if (is_view_)
                {
                    data_ = other.data_;
                }
                else
                {
                    data_ = new T[size_];
                    std::copy(other.data_, other.data_+size_, data_);
                }
            }

            const_marray(const_marray&& other)
            : data_(other.data_), size_(other.size_), stride_(other.stride_),
              is_view_(other.is_view_)
            {
                if (!is_view_)
                {
                    other.data_ = NULL;
                    other.reset();
                }
            }

            const_marray(const const_marray& other, const construct_view_t& cv)
            : data_(other.data_), size_(other.size_), stride_(other.stride_),
              is_view_(true) {}

            explicit const_marray(idx_type n = 0, const T& val = T())
            : data_(NULL), size_(n), stride_(1), is_view_(false)
            {
                if (size_ > 0)
                {
                    data_ = new T[size_];
                    std::fill(data_, data_+size_, val);
                }
            }

            template <unsigned ndim_>
            const_marray(const const_marray_ref<T, ndim_, ndim_>& other)
            : data_(other.array_.data_+other.idx), size_(other.array_.len_[ndim_-1]),
              stride_(other.array_.stride_[ndim_-1]), is_view_(true) {}

            ~const_marray()
            {
                if (!is_view_ && data_) delete[] data_;
            }

            bool isView() const
            {
                return is_view_;
            }

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
                assert(size_ > 0);
                return *data_;
            }

            const_reference back() const
            {
                assert(size_ > 0);
                return *(data_+(size_-1)*stride_);
            }

            const_marray& operator=(const const_marray& other) = delete;

            const_reference operator[](idx_type i) const
            {
                assert(i >= 0 && i < size_);
                return *(data_+i*stride_);
            }

            const_marray<T, 1>& operator[](const slice::slice_all_t& x) const
            {
                return *this;
            }

            const_reference operator()(idx_type i) const
            {
                return (*this)[i];
            }

            const_marray<T, 1>& operator()(const slice::slice_all_t& x) const
            {
                return *this;
            }

            const_pointer data() const
            {
                return data_;
            }

            idx_type length() const
            {
                return size_;
            }

            idx_type length(unsigned dim) const
            {
                return size_;
            }

            std::vector<idx_type> lengths() const
            {
                return vec(size_);
            }

            size_type size() const
            {
                return size_;
            }
    };

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
            using const_marray<T, 1>::stride_;
            using const_marray<T, 1>::is_view_;

            marray(const const_marray<T, 1>& other)
            : const_marray<T, 1>(other) {}

        public:
            marray(marray& other)
            : const_marray<T, 1>(other) {}

            marray(marray&& other)
            : const_marray<T, 1>(move(other)) {}

            marray(marray& other, const construct_view_t& cv)
            : const_marray<T, 1>(other, cv) {}

            explicit marray(idx_type n = 0, const T& val = T())
            :  const_marray<T, 1>(n, val) {}

            template <unsigned ndim_>
            marray(const const_marray_ref<T, ndim_, ndim_>& other)
            : const_marray<T, 1>(other) {}

            using const_marray<T, 1>::reset;

            void clear()
            {
                reset();
            }

            void unView()
            {
                if (!is_view_) return;

                const_marray<T, 1> old(*this);
                reset(size_);
                copy(old, *this);
            }

            void push_back(unsigned dim, const T& x)
            {
                assert(dim == 0);
                push_back(x);
            }

            void pop_back(unsigned dim)
            {
                assert(dim == 0);
                pop_back();
            }

            void push_back(const T& x)
            {
                unView();
                pointer old_data = data_;
                data_ = new T[++size_];
                std::copy(old_data, old_data+size_-1, data_);
                data_[size_-1] = x;
            }

            void pop_back()
            {
                assert(size_ > 0);
                unView();
                size_--;
            }

            reference front(unsigned dim)
            {
                return const_cast<reference>(const_marray<T, 1>::front(dim));
            }

            reference back(unsigned dim)
            {
                return const_cast<reference>(const_marray<T, 1>::back(dim));
            }

            reference front()
            {
                return const_cast<reference>(const_marray<T, 1>::front());
            }

            reference back()
            {
                return const_cast<reference>(const_marray<T, 1>::back());
            }

            marray& operator=(const marray& other)
            {
                copy(other, *this);
                return *this;
            }

            marray& operator=(const const_marray<T, 1>& other)
            {
                copy(other, *this);
                return *this;
            }

            marray& operator=(const T& x)
            {
                for (size_type i = 0;i < size_;i++)
                {
                    *(data_+i*stride_) = x;
                }
                return *this;
            }

            using const_marray<T, 1>::operator[];

            reference operator[](idx_type i)
            {
                assert(i >= 0 && i < size_);
                return *(data_+i*stride_);
            }

            marray<T, 1>& operator[](const slice::slice_all_t& x)
            {
                return *this;
            }

            using const_marray<T, 1>::operator();

            reference operator()(idx_type i)
            {
                return (*this)[i];
            }

            marray<T, 1>& operator()(const slice::slice_all_t& x)
            {
                return *this;
            }

            pointer data()
            {
                return data_;
            }

            void swap(marray&& other)
            {
                swap(other);
            }

            void swap(marray& other)
            {
                using std::swap;
                swap(data_,      other.data_);
                swap(size_,      other.size_);
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

            friend void copy(const const_marray<T, 1>& a, marray&& b)
            {
                copy(a, b);
            }

            friend void copy(const const_marray<T, 1>& a, marray& b)
            {
                assert(a.size_ == b.size_);
                const_pointer a_ = a.data();
                      pointer b_ = b.data();
                for (unsigned i = 0;i < a.size_;i++)
                {
                    *b_ = *a_;
                    a_ += a.stride_;
                    b_ += b.stride_;
                }
            }
    };

    template <typename T, unsigned ndim, unsigned dim>
    void copy(const const_marray<T, ndim-dim-1>& a, marray_ref<T, ndim, dim>&& b)
    {
        marray<T, ndim-dim-1> b_(b);
        copy(a, b_);
    }

    template <typename T, unsigned ndim, unsigned newdim, unsigned dim>
    void copy(const const_marray<T, ndim+newdim-dim-1>& a, marray_slice<T, ndim, newdim, dim>&& b)
    {
        marray<T, ndim+newdim-dim-1> b_(b);
        copy(a, b_);
    }

    template <typename T>
    class const_row : public const_marray<T, 1>
    {
        public:
            using typename const_marray<T, 1>::value_type;
            using typename const_marray<T, 1>::idx_type;
            using typename const_marray<T, 1>::size_type;
            using typename const_marray<T, 1>::reference;
            using typename const_marray<T, 1>::const_reference;
            using typename const_marray<T, 1>::pointer;
            using typename const_marray<T, 1>::const_pointer;

            const_row() {}

            const_row(const const_row& other) : const_marray<T, 1>(other) {}

            const_row(const_row&& other) : const_marray<T, 1>(move(other)) {}

            const_row(const const_row& other, const construct_view_t& cv) : const_marray<T, 1>(other, cv) {}

            const_row(const const_marray<T, 1>& other) : const_marray<T, 1>(other) {}

            const_row(const_marray<T, 1>&& other) : const_marray<T, 1>(move(other)) {}

            const_row(const const_marray<T, 1>& other, const construct_view_t& cv) : const_marray<T, 1>(other, cv) {}

            explicit const_row(idx_type m, const T& val = T())
            : const_marray<T, 1>(m, val) {}

            const_row& operator=(const const_row& other) = delete;
    };

    template <typename T>
    class row : public marray<T, 1>
    {
        public:
            using typename marray<T, 1>::value_type;
            using typename marray<T, 1>::idx_type;
            using typename marray<T, 1>::size_type;
            using typename marray<T, 1>::reference;
            using typename marray<T, 1>::const_reference;
            using typename marray<T, 1>::pointer;
            using typename marray<T, 1>::const_pointer;

            using marray<T, 1>::operator=;

            row() {}

            row(row& other) : marray<T, 1>(other) {}

            row(row&& other) : marray<T, 1>(move(other)) {}

            row(row& other, const construct_view_t& cv) : marray<T, 1>(other, cv) {}

            row(marray<T, 1>& other) : marray<T, 1>(other) {}

            row(marray<T, 1>&& other) : marray<T, 1>(move(other)) {}

            row(marray<T, 1>& other, const construct_view_t& cv) : marray<T, 1>(other, cv) {}

            explicit row(idx_type m, const T& val = T())
            : marray<T, 1>(m, val) {}

            row& operator=(const row& other)
            {
                copy(other, *this);
                return *this;
            }
    };

    template <typename T>
    class const_matrix : public const_marray<T, 2>
    {
        public:
            using typename const_marray<T, 2>::value_type;
            using typename const_marray<T, 2>::idx_type;
            using typename const_marray<T, 2>::size_type;
            using typename const_marray<T, 2>::reference;
            using typename const_marray<T, 2>::const_reference;
            using typename const_marray<T, 2>::pointer;
            using typename const_marray<T, 2>::const_pointer;

            const_matrix() {}

            const_matrix(const const_matrix& other) : const_marray<T, 2>(other) {}

            const_matrix(const_matrix&& other) : const_marray<T, 2>(move(other)) {}

            const_matrix(const const_matrix& other, const construct_view_t& cv) : const_marray<T, 2>(other, cv) {}

            const_matrix(const const_marray<T, 2>& other) : const_marray<T, 2>(other) {}

            const_matrix(const_marray<T, 2>&& other) : const_marray<T, 2>(move(other)) {}

            const_matrix(const const_marray<T, 2>& other, const construct_view_t& cv) : const_marray<T, 2>(other, cv) {}

            const_matrix(idx_type m, idx_type n, const T& val = T(), Layout layout = DEFAULT)
            : const_marray<T, 2>({m,n}, val, layout) {}

            const_matrix& operator=(const const_matrix& other) = delete;
    };

    template <typename T>
    class matrix : public marray<T, 2>
    {
        public:
            using typename marray<T, 2>::value_type;
            using typename marray<T, 2>::idx_type;
            using typename marray<T, 2>::size_type;
            using typename marray<T, 2>::reference;
            using typename marray<T, 2>::const_reference;
            using typename marray<T, 2>::pointer;
            using typename marray<T, 2>::const_pointer;

            using marray<T, 2>::operator=;

            matrix() {}

            matrix(matrix& other) : marray<T, 2>(other) {}

            matrix(matrix&& other) : marray<T, 2>(move(other)) {}

            matrix(matrix& other, const construct_view_t& cv) : marray<T, 2>(other, cv) {}

            matrix(marray<T, 2>& other) : marray<T, 2>(other) {}

            matrix(marray<T, 2>&& other) : marray<T, 2>(move(other)) {}

            matrix(marray<T, 2>& other, const construct_view_t& cv) : marray<T, 2>(other, cv) {}

            matrix(idx_type m, idx_type n, const T& val = T(), Layout layout = DEFAULT)
            : marray<T, 2>({m,n}, val, layout) {}

            void reset(idx_type m, idx_type n, const T& val = T(), Layout layout = DEFAULT)
            {
                marray<T, 2>::reset({m,n}, val, layout);
            }

            matrix& operator=(const matrix& other)
            {
                copy(other, *this);
                return *this;
            }
    };

}

#endif
