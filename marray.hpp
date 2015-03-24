#ifndef _MARRAY_MARRAY_HPP_
#define _MARRAY_MARRAY_HPP_

#ifndef MARRAY_TEST
#define MARRAY_TEST(...)
#endif

#include <type_traits>
#include <array>
#include <cassert>
#include <algorithm>

#include "iterator.hpp"
#include "utility.hpp"

namespace MArray
{
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
                               const std::array<unsigned,newdim-1>& dims,
                               const std::array<idx_type,newdim-1>& lens, idx_type i)
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
                         const std::array<unsigned,newdim-1>& dims,
                         const std::array<idx_type,newdim-1>& lens, idx_type i)
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
                assert(i >= 0 && i < array_.len[dim-1]);
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
     * into an array view (of dimension ndim-dim+1+newdim) or indexed one last
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
                               const std::array<unsigned,newdim-1>& dims,
                               const std::array<idx_type,newdim-1>& lens, idx_type i)
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
                         const std::array<unsigned,newdim-1>& dims,
                         const std::array<idx_type,newdim-1>& lens, idx_type i)
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

            template <typename I>
            marray<T, newdim+1> operator[](const range_t<I>& x)
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= array_.len_[ndim-1]);
                marray<T, newdim+1> ret;
                ret.data_ = array_.data_+idx+array_.stride[ndim-1]*x.front();
                ret.size_ = 1;
                ret.is_view_ = true;
                ret.layout_ = array_.layout_;
                for (unsigned dim = 0;dim < newdim;dim++)
                {
                    ret.len_[dim] = array_.len_[dims[dim]];
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
     * Specifies the layout of the array data. The layout is only significant
     * when a data
     */
    enum Layout {COLUMN_MAJOR, ROW_MAJOR, DEFAULT=ROW_MAJOR};

    namespace detail
    {
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

        /*
         * These helper structs determine the return type after indexing an
         * array with the given index types. The return type may be either
         * a reference to the underlying datatype (no slicing), or a lower-
         * dimensional view of the array if one or more arguments are ranges.
         */
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


        /*
         * These helper structs are identical to const_return_type except
         * that a non-const type is returned.
         */
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

        /*
         * These helper structs apply the given set of indices to an array
         * to generate either a reference to a specific element or a
         * lower-dimensional view of the array.
         */
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


        /*
         * These helper structs are identical to get_const_slice except that a
         * non-const reference is returned.
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

        /*
         * Helper class to determine if the arguments Args... are of one of the
         * following forms which are acceptable for constructing an instance of
         * [const_]marray<T, N>:
         *
         * 1) Args... = <int1, int2, ..., intN>
         *
         * 2) Args... = <int1, int2, ..., intN, T>
         *
         * 3) Args... = <int1, int2, ..., intN, T, Layout>
         *
         * 4) Args... = <int1, int2, ..., intN, [const] T*>
         *
         * 5) Args... = <int1, int2, ..., intN, [const] T*, Layout>
         *
         * where intX is any integral type. The additional forms allow for default
         * values of the fill element and layout parameters (1-3) and for wrapping
         * an external pointer (4-5).
         *
         * The default definition matches any cases which do not conform to one
         * of the above forms and derives from false_type.
         */
        template <typename T, unsigned ndim, unsigned dim, typename... Args>
        struct are_marray_args : std::false_type {};

        /*
         * This specialization checks whether the argument at position I-1 is
         * integral or not. If so, it recursively derives from are_marray_args
         * to check the remaining arguments. If not, it derives from false_type.
         */
        template <typename T, unsigned ndim, unsigned dim, typename U, typename... Args>
        struct are_marray_args<T, ndim, dim, U, Args...>
        : std::conditional<std::is_integral<typename std::remove_reference<U>::type>::value,
                           are_marray_args<T, ndim, dim+1, Args...>,
                           std::false_type>::type {};

        /*
         * This specialization is reached when the first I arguments are integral
         * and there are no remaining arguments (type 1 above).
         */
        template <typename T, unsigned ndim>
        struct are_marray_args<T, ndim, ndim> : std::true_type {};

        /*
         * This specialization is reached when the first I arguments are integral
         * and the only remaining argument is convertible to T (type 2 above) or
         * to const T* (type 4).
         */
        template <typename T, unsigned ndim, typename U>
        struct are_marray_args<T, ndim, ndim, U>
        : std::integral_constant<bool,std::is_convertible<U,T>::value ||
                                      std::is_convertible<U,const T*>::value> {};

        /*
         * This specialization is reached when the first I arguments are integral
         * and the two remaining arguments are convertible to T and of type
         * Layout respectively (type 3 above) or convertible to const T* and of
         * type Layout (type 5).
         */
        template <typename T, unsigned ndim, typename U>
        struct are_marray_args<T, ndim, ndim, U, Layout>
        : std::integral_constant<bool,std::is_convertible<U,T>::value ||
                                      std::is_convertible<U,const T*>::value> {};

        /*
         * This helper class initializes a const_marray<T, ndim> from one of the
         * argument forms allowed by are_marray_args above. The default definition
         * consumes one integral argument and assigns it to the (dim-1)th length
         * parameter of the array.
         */
        template <typename T, unsigned ndim, unsigned dim, typename U, typename... Args>
        struct marray_reset
        {
            marray_reset(const_marray<T, ndim>& array, U len, Args&&... args)
            {
                array.len_[dim-1] = len;
                marray_reset<T, ndim, dim+1, Args...>(array, std::forward<Args>(args)...);
            }
        };

        /*
         * This specialization is triggered when all integral arguments have been
         * consumed. Any remaining arguments are passed to const_marray::reset
         * to complete the initialization.
         */
        template <typename T, unsigned ndim, typename U, typename... Args>
        struct marray_reset<T, ndim, ndim, U, Args...>
        {
            marray_reset(const_marray<T, ndim>& array, U len, Args&&... args)
            {
                array.len_[ndim-1] = len;
                array.reset(array.len_, std::forward<Args>(args)...);
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
     * 2) explicit const_marray<T, ndim>(const std::array<idx_type, ndim>& len,
     *                                  const T& val = T(), Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout. The array
     *  is initialized with the element val.
     *
     * 3) template <typename I>
     *    explicit const_marray<T, ndim>(std::initializer_list<I> len,
     *                                   const T& val = T(), Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout. The array
     *  is initialized with the element val. The initializer list must have
     *  exactly ndim entries.
     *
     * 4) template <typename I1, ..., typename IN>
     *    explicit const_marray<T, ndim>(I1 len1, ..., IN lenN,
     *                                   const T& val = T(), Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout. The array
     *  is initialized with the element val. There must be exactly ndim lengths
     *  specified, and each must have a (possibly different) integral type.
     *
     * 5) const_marray<T, ndim>(const std::array<idx_type, ndim>& len,
     *                          const T* ptr, Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout, using ptr
     *  as the array data. The given pointer is referenced directly by the
     *  array, providing a "view" of an externally-allocated data segment.
     *  No internal allocation or copying is performed.
     *
     * 6) template <typename I>
     *    const_marray<T, ndim>(std::initializer_list<I> len,
     *                          const T* ptr, Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout, using ptr
     *  as the array data. The given pointer is referenced directly by the
     *  array, providing a "view" of an externally-allocated data segment.
     *  No internal allocation or copying is performed. The initializer list
     *  must have exactly ndim entries.
     *
     * 7) template <typename I>
     *    const_marray<T, ndim>(I len1, ..., I lenN,
     *                          const T* ptr, Layout layout = DEFAULT)
     *
     *  Constructs an array with the given lengths and data layout, using ptr
     *  as the array data. The given pointer is referenced directly by the
     *  array, providing a "view" of an externally-allocated data segment.
     *  No internal allocation or copying is performed. There must be exactly
     *  ndim lengths specified, and each must have a (possibly different)
     *  integral type.
     *
     * 8) const_marray<T, ndim>(const const_marray<T, ndim>& other)
     *
     *  Creates a copy of the given array. If the array is a view, then so is
     *  the copy.
     *
     * 9) const_marray<T, ndim>(const_marray<T, ndim>&& other)
     *
     *  Move-constructs a copy of the given array. If the array is a view, then
     *  this behaves identically to the copy constructor. Otherwise, the given
     *  array is reinitialized to a state as if it were constructed using the
     *  empty constructor.
     *
     * 10) const_marray<T, ndim>(const const_marray<T, ndim>& other,
     *                           const construct_view_t& cv)
     *
     *  Constructs a view of the given array, even if it is not a view itself.
     *
     * 11) const_marray<T, ndim>(const const_marray<T, ndim>& other,
     *                           const construct_copy_t& cp)
     *
     *  Constructs a new (freshly allocated) copy of the given array, even if
     *  it is a view.
     *
     * 12) template <unsigned dim>
     *     const_marray<T, ndim>(const const_marray_ref<T, ndim+dim-1, dim>& ref)
     *
     *  Creates a view from the given array reference as if it were repeatedly
     *  indexed with slice::all.
     *
     * 13) template <unsigned dim, unsigned newdim>
     *     const_marray<T, ndim>(const const_marray_slice<T, ndim+newdim+dim-1, dim>& slice)
     *
     *  Creates a view from the given array slice as if it were repeatedly
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
        template <typename T_, unsigned ndim_, unsigned dim_, typename U, typename... Args> friend class detail::marray_reset;

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

            void reset()
            {
                if (!is_view_ && data_) delete[] data_;
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
                assert(a.size() == 0);
                assert(a.lengths() == make_array<idx_type>(0, 0, 0));
                assert(a.strides() == make_array<size_type>(0, 0, 0));
                assert(a.layout() == DEFAULT);
            )

            void reset(const std::array<idx_type, ndim>& len,
                       const T& val = T(), Layout layout = DEFAULT)
            {
                size_type old_size = (is_view_ ? 0 : size_);
                len_ = len;

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

                if (size_ > old_size || is_view_ || !data_)
                {
                    if (!is_view_ && data_) delete[] data_;
                    data_ = new T[size_];
                }
                std::fill(data_, data_+size_, val);

                layout_ = layout;
                is_view_ = false;
            }

            MARRAY_TEST
            (
                marray<double, 3> a;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                a.reset(make_array<idx_type>(2, 4, 5));
                assert(a.data() != NULL);
                assert(a.size() == 2*4*5);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(4*5, 5, 1));
                assert(a.layout() == DEFAULT);
                a.reset();
            )

            void reset(const std::array<idx_type, ndim>& len,
                       pointer ptr, Layout layout = DEFAULT)
            {
                size_type old_size = (is_view_ ? 0 : size_);
                len_ = len;

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

                if (!is_view_ && data_) delete[] data_;
                data_ = ptr;

                layout_ = layout;
                is_view_ = true;
            }

            MARRAY_TEST
            (
                double p;
                marray<double, 3> a;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                a.reset(make_array<idx_type>(2, 4, 5), &p);
                assert(a.data() == &p);
                assert(a.size() == 2*4*5);
                assert(a.lengths() == make_array<idx_type>(2, 4, 5));
                assert(a.strides() == make_array<size_type>(4*5, 5, 1));
                assert(a.layout() == DEFAULT);
                // make sure this doesn't segfault
                a.reset();
            )

            template <typename... Args>
            void reset(typename std::enable_if<detail::are_marray_args<T, ndim, 1, Args...>::value,idx_type>::type len0, Args&&... args)
            {
                detail::marray_reset<T, ndim, 1, idx_type, Args...>(*this, len0, std::forward<Args>(args)...);
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

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 3> ca(a);
                const_marray<double, 3> va(a, construct_view);
                const_marray<double, 3> cva(va);

                assert(ca.data() != NULL);
                assert(ca.data() != a.data());
                assert(ca.size() == a.size());
                assert(ca.lengths() == a.lengths());
                assert(ca.strides() == a.strides());
                assert(ca.layout() == a.layout());
                assert(!ca.isView());

                assert(cva.data() == a.data());
                assert(cva.size() == a.size());
                assert(cva.lengths() == a.lengths());
                assert(cva.strides() == a.strides());
                assert(cva.layout() == a.layout());
                assert(cva.isView());
            )

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

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                const double *data = a.data();
                const_marray<double, 3> ca(std::move(a));

                assert(ca.data() == data);
                assert(ca.size() == 2*4*5);
                assert(ca.lengths() == make_array<idx_type>(2, 4, 5));
                assert(ca.strides() == make_array<size_type>(1, 2, 2*4));
                assert(ca.layout() == COLUMN_MAJOR);
                assert(!ca.isView());

                assert(a.data() == NULL);
                assert(a.size() == 0);
                assert(a.lengths() == make_array<idx_type>(0, 0, 0));
                assert(a.strides() == make_array<size_type>(0, 0, 0));
                assert(a.layout() == DEFAULT);
            )

            const_marray(const const_marray& other, const construct_view_t& cv)
            : data_(other.data_), size_(other.size_), len_(other.len_),
              stride_(other.stride_), is_view_(true),
              layout_(other.layout_) {}

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 3> va(a, construct_view);

                assert(va.data() == a.data());
                assert(va.size() == a.size());
                assert(va.lengths() == a.lengths());
                assert(va.strides() == a.strides());
                assert(va.layout() == a.layout());
                assert(va.isView());
            )

            const_marray(const const_marray& other, const construct_copy_t& cp)
            : data_(NULL), size_(other.size_), len_(other.len_),
              stride_(other.stride_), is_view_(false),
              layout_(other.layout_)
            {
                if (size_ > 0)
                {
                    data_ = new T[size_];
                    copy(other, static_cast<marray<T, ndim>&>(*this));
                }
            }

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5, 0.0, COLUMN_MAJOR);
                const_marray<double, 3> va(a, construct_view);
                const_marray<double, 3> ca(va, construct_copy);

                assert(ca.data() != NULL);
                assert(ca.data() != a.data());
                assert(ca.size() == a.size());
                assert(ca.lengths() == a.lengths());
                assert(ca.strides() == a.strides());
                assert(ca.layout() == a.layout());
                assert(!ca.isView());
            )

            const_marray()
            : data_(NULL), size_(0), len_(), stride_(), is_view_(false), layout_(DEFAULT) {}

            explicit const_marray(const std::array<idx_type, ndim>& len,
                                  const T& val = T(), Layout layout = DEFAULT)
            : data_(NULL), size_(0), len_(len), stride_(), is_view_(false), layout_(layout)
            {
                reset(len_, val, layout_);
            }

            template <typename U>
            explicit const_marray(std::initializer_list<U> len, const T& val = T(),
                                  Layout layout = DEFAULT, typename std::enable_if<std::is_integral<U>::value>::type* = 0)
            : data_(NULL), size_(0), len_(), stride_(), is_view_(false), layout_(layout)
            {
                assert(len.size() == ndim);
                std::copy_n(len.begin(), ndim, len_.begin());
                reset(len_, val, layout_);
            }

            MARRAY_TEST
            (
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                const_marray<double, 3> a(make_array<idx_type>(2, 4, 5));
                const_marray<double, 3> b({2, 4, 5});

                assert(a.data() != NULL);
                assert(b.data() != NULL);
                assert(a.size() == b.size());
                assert(a.lengths() == b.lengths());
                assert(a.strides() == b.strides());
            )

            explicit const_marray(const std::array<idx_type, ndim>& len,
                                  const_pointer ptr, Layout layout = DEFAULT)
            : data_(NULL), size_(0), len_(len), stride_(), is_view_(false), layout_(layout)
            {
                reset(len_, const_cast<pointer>(ptr), layout_);
            }

            template <typename U>
            const_marray(std::initializer_list<U> len, const_pointer ptr,
                         Layout layout = DEFAULT, typename std::enable_if<std::is_integral<U>::value>::type* = 0)
            : data_(NULL), size_(0), len_(), stride_(), is_view_(false), layout_(layout)
            {
                assert(len.size() == ndim);
                std::copy_n(len.begin(), ndim, len_.begin());
                reset(len_, const_cast<pointer>(ptr), layout_);
            }

            MARRAY_TEST
            (
                double p;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                const_marray<double, 3> a(make_array<idx_type>(2, 4, 5), &p);
                const_marray<double, 3> b({2, 4, 5}, &p);

                assert(a.data() == &p);
                assert(b.data() == &p);
                assert(a.size() == b.size());
                assert(a.lengths() == b.lengths());
                assert(a.strides() == b.strides());
            )

            template <typename... Args>
            explicit const_marray(typename std::enable_if<detail::are_marray_args<T, ndim, 1, Args...>::value,idx_type>::type len0, Args&&... args)
            : data_(NULL), size_(0), is_view_(false), layout_(DEFAULT)
            {
                detail::marray_reset<T, ndim, 1, idx_type, Args...>(*this, len0, std::forward<Args>(args)...);
            }

            MARRAY_TEST
            (
                const_marray<double, 3> a(2, 4, 5);
                const_marray<double, 3> b(2, 4, 5, 0.0, DEFAULT);

                assert(a.data() != NULL);
                assert(b.data() != NULL);
                assert(a.size() == b.size());
                assert(a.lengths() == b.lengths());
                assert(a.strides() == b.strides());
            )

            MARRAY_TEST
            (
                double p;
                const_marray<double, 3> a(2, 4, 5, &p);
                const_marray<double, 3> b(2, 4, 5, &p, DEFAULT);

                assert(a.data() == &p);
                assert(b.data() == &p);
                assert(a.size() == b.size());
                assert(a.lengths() == b.lengths());
                assert(a.strides() == b.strides());
            )

            template <unsigned ndim_>
            const_marray(const const_marray_ref<T, ndim_, ndim_-ndim+1>& other)
            : data_(other.array_.data_+other.idx),
              is_view_(true), layout_(other.array_.layout_)
            {
                std::copy(other.array_.len_.begin()+ndim_-ndim, other.array_.len_.end(), len_.begin());
                std::copy(other.array_.stride_.begin()+ndim_-ndim, other.array_.stride_.end(), stride_.begin());

                size_ = 1;
                for (idx_type i = 0;i < ndim;i++) size_ *= len_[i];
            }

            MARRAY_TEST
            (
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                const_marray<double, 3> a(2, 4, 5, COLUMN_MAJOR);
                const_marray<double, 2> b = a[0];

                assert(a.data() == b.data());
                assert(b.size() == 4*5);
                assert(b.lengths() == make_array<idx_type>(4, 5));
                assert(b.strides() == make_array<size_type>(2, 2*4));
                assert(b.layout() == COLUMN_MAJOR);
            )

            template <unsigned ndim_, unsigned newdim_>
            const_marray(const const_marray_slice<T, ndim_, ndim_-ndim+1+newdim_, newdim_>& other)
            : data_(other.array_.data_+other.idx),
              is_view_(true), layout_(other.array_.layout_)
            {
                for (unsigned i = 0;i < newdim_;i++)
                {
                    len_[i] = other.lens[i];
                    stride_[i] = other.array_.stride_[other.dims[i]];
                }
                std::copy(other.array_.len_.begin()+ndim_-ndim+newdim_, other.array_.len_.end(), len_.begin()+newdim_);
                std::copy(other.array_.stride_.begin()+ndim_-ndim+newdim_, other.array_.stride_.end(), stride_.begin()+newdim_);

                size_ = 1;
                for (idx_type i = 0;i < ndim;i++) size_ *= len_[i];
            }

            MARRAY_TEST
            (
                using slice::all;
                typedef marray<double, 3>::idx_type idx_type;
                typedef marray<double, 3>::size_type size_type;
                const_marray<double, 4> a(2, 1, 4, 5, COLUMN_MAJOR);
                const_marray<double, 3> b = a[0][all];

                assert(a.data() == b.data());
                assert(b.size() == 1*4*5);
                assert(b.lengths() == make_array<idx_type>(1, 4, 5));
                assert(b.strides() == make_array<size_type>(2, 2*2, 2*1*4));
                assert(b.layout() == COLUMN_MAJOR);
            )

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

            template <typename I>
            const_marray_slice<T, ndim, 2, 1> operator[](const range_t<I>& x) const
            {
                assert(x.fron() >= 0 && x.back() <= len_[0]);
                return const_marray_slice<T, ndim, 2, 1>(*this, (size_type)0, {}, {}, x);
            }

            const_marray_slice<T, ndim, 2, 1> operator[](const slice::all_t& x) const
            {
                return const_marray_slice<T, ndim, 2, 1>(*this, (size_type)0, {}, {}, range(idx_type(), len_[0]));
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

            size_type size() const
            {
                return size_;
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
            : const_marray<T, ndim>(std::move(other)) {}

            marray(marray& other, const construct_view_t& cv)
            : const_marray<T, ndim>(other, cv) {}

            marray(const marray& other, const construct_copy_t& cp)
            : const_marray<T, ndim>(other, cp) {}

            marray() {}

            explicit marray(const std::array<idx_type, ndim>& len, const T& val = T(), Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, val, layout) {}

            explicit marray(std::initializer_list<idx_type> len, const T& val = T(), Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, val, layout) {}

            explicit marray(const std::array<idx_type, ndim>& len, pointer ptr, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, ptr, layout) {}

            explicit marray(std::initializer_list<idx_type> len, pointer ptr, Layout layout = DEFAULT)
            : const_marray<T, ndim>(len, ptr, layout) {}

            template <typename... Args>
            marray(typename std::enable_if<detail::are_marray_args<T, ndim, 1, Args...>::value,idx_type>::type len0, Args&&... args)
            : const_marray<T, ndim>(len0, std::forward<Args>(args)...) {}

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

            void resize(const std::array<idx_type, ndim>& len, const T& val = T(), Layout layout = DEFAULT)
            {
                if (layout == DEFAULT) layout = layout_;

                std::array<idx_type, ndim> common;
                for (unsigned i = 0;i < ndim;i++) common[i] = std::min(len_[i], len[i]);

                const_marray<T, ndim> old(std::move(*this));
                reset(len, val, layout);

                if (old.data_ && data_)
                {
                    Iterator<idx_type, size_type> it(common, old.stride_, stride_);
                    const_pointer a_ = old.data_;
                          pointer b_ =     data_;
                    while (it.nextIteration(a_, b_)) *b_ = *a_;
                }
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

                std::array<size_type, ndim-1> slice_stride_;
                std::copy_n(stride_.begin(), dim, slice_stride_.begin());
                std::copy_n(stride_.begin()+dim+1, ndim-dim, slice_stride_.begin()+dim);

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
                *this = static_cast<const const_marray<T, ndim>&>(other);
                return *this;
            }

            marray& operator=(const const_marray<T, ndim>& other)
            {
                if (other.len_ != len_) reset(other.len_, T(), other.layout_);

                if (!other.data_ || !data_) return *this;

                if (!other.is_view_ && !is_view_ && other.layout_ == layout_)
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
                assert(a.lengths() == b.lengths());
                b = a;
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
            pointer data_;
            idx_type size_;
            size_type stride_;
            bool is_view_;

            void reset(idx_type n = 0, const T& val = T())
            {
                idx_type old_size = (is_view_ ? 0 : size_);
                size_ = n;

                if (size_ > old_size || is_view_ || !data_)
                {
                    if (!is_view_ && data_) delete[] data_;
                    data_ = new T[size_];
                }
                std::fill(data_, data_+size_, val);

                is_view_ = false;
            }

            void reset(idx_type n, pointer ptr)
            {
                idx_type old_size = (is_view_ ? 0 : size_);
                size_ = n;

                if (!is_view_ && data_) delete[] data_;
                data_ = ptr;

                is_view_ = true;
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

            const_marray(const const_marray& other, const construct_copy_t& cp)
            : data_(NULL), size_(other.size_), stride_(other.stride_),
              is_view_(false)
            {
                data_ = new T[size_];
                copy(other, static_cast<marray<T, 1>&>(*this));
            }

            explicit const_marray(idx_type n = 0, const T& val = T())
            : data_(NULL), size_(n), stride_(1), is_view_(false)
            {
                reset(size_, val);
            }

            explicit const_marray(idx_type n, const_pointer ptr)
            : data_(const_cast<pointer>(ptr)), size_(n), stride_(1), is_view_(true) {}

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

            template <typename I>
            const_marray<T, 1> operator[](const range_t<I>& x) const
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= size_);
                const_marray<T, 1> ret(*this, construct_view);
                ret.size_ = x.size();
                ret.data_ += stride_*x.front();
                return ret;
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

            std::array<idx_type, 1> lengths() const
            {
                return make_array(size_);
            }

            size_type stride() const
            {
                return stride_;
            }

            size_type stride(unsigned dim) const
            {
                return stride_;
            }

            std::array<size_type, 1> strides() const
            {
                return make_array(stride_);
            }

            size_type size() const
            {
                return size_;
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

            marray(const marray& other, const construct_copy_t& cp)
            : const_marray<T, 1>(other, cp) {}

            explicit marray(idx_type n = 0, const T& val = T())
            : const_marray<T, 1>(n, val) {}

            explicit marray(idx_type n, const_pointer ptr)
            : const_marray<T, 1>(n, ptr) {}

            template <unsigned ndim_>
            marray(const const_marray_ref<T, ndim_, ndim_>& other)
            : const_marray<T, 1>(other) {}

            using const_marray<T, 1>::reset;

            void clear()
            {
                reset();
            }

            void resize(idx_type n, const T& val = T())
            {
                size_type common = std::min(n, size_);

                const_marray<T, 1> old(std::move(*this));
                reset(n, val);

                const_pointer a_ = old.data_;
                      pointer b_ =     data_;
                for (size_type i = 0;i < common;i++)
                {
                    *b_ = *a_;
                    a_ += old.stride_;
                    b_ +=     stride_;
                }
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

            template <typename I>
            marray<T, 1> operator[](const range_t<I>& x)
            {
                assert(x.front() <= x.back() && x.front() >= 0 && x.back() <= size_);
                marray<T, 1> ret(*this, construct_view);
                ret.size_ = x.size();
                ret.data_ += stride_*x.front();
                return ret;
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
                assert(a.size() == b.size());
                const_pointer a_ = a.data();
                      pointer b_ = b.data();
                for (size_type i = 0;i < a.size();i++)
                {
                    *b_ = *a_;
                    a_ += a.stride();
                    b_ += b.stride();
                }
            }
    };

    /*
     * These helper functions are required to facilitate copying into temporary
     * views created by incompletely indexing an array (since implicit
     * conversion will not occur for rvalue references).
     */
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

    /*
     * Convenient names for 1- and 2-dimensional array types.
     */
    template <typename T> using const_row = const_marray<T, 1>;
    template <typename T> using row = marray<T, 1>;
    template <typename T> using const_matrix = const_marray<T, 2>;
    template <typename T> using matrix = marray<T, 2>;
}

#endif
