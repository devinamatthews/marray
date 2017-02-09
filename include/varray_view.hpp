#ifndef _MARRAY_VARRAY_VIEW_HPP_
#define _MARRAY_VARRAY_VIEW_HPP_

#include "varray.hpp"

namespace MArray
{

template <typename T>
class varray_view
{
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;

    protected:
        pointer data_ = nullptr;
        std::vector<idx_type> len_;
        std::vector<stride_type> stride_;

        template <unsigned Dim>
        void get_slice(pointer& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride) const {}

        template <unsigned Dim, typename... Args>
        void get_slice(pointer& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride,
                       idx_type arg, Args&&... args) const
        {
            MARRAY_ASSERT(arg >= 0 && arg < len_[Dim]);
            ptr += arg*stride_[Dim];
            get_slice<Dim+1>(ptr, len, stride, std::forward<Args>(args)...);
        }

        template <unsigned Dim, typename I, typename... Args>
        void get_slice(pointer& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride,
                       const range_t<I>& arg, Args&&... args) const
        {
            MARRAY_ASSERT(arg.front() <= arg.back());
            MARRAY_ASSERT(arg.front() >= 0 && arg.back() < len_[Dim]);
            ptr += arg.front()*stride_[Dim];
            len.push_back(arg.size());
            stride.push_back(arg.step()*stride_[Dim]);
            get_slice<Dim+1>(ptr, len, stride, std::forward<Args>(args)...);
        }

        template <unsigned Dim, typename... Args>
        void get_slice(pointer& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride,
                       all_t, Args&&... args) const
        {
            len.push_back(len_[Dim]);
            stride.push_back(stride_[Dim]);
            get_slice<Dim+1>(ptr, len, stride, std::forward<Args>(args)...);
        }

        template <unsigned Dim>
        void get_reference(pointer& ptr) const {}

        template <unsigned Dim, typename... Args>
        void get_reference(pointer& ptr, idx_type arg, Args&&... args) const
        {
            MARRAY_ASSERT(arg >= 0 && arg < len_[Dim]);
            ptr += arg*stride_[Dim];
            get_reference<Dim+1>(ptr, std::forward<Args>(args)...);
        }

    public:
        varray_view() {}

        varray_view(const varray_view& other)
        {
            reset(other);
        }

        varray_view(varray_view&& other)
        {
            reset(std::move(other));
        }

        template <typename U, typename=detail::enable_if_convertible_t<U*,pointer>>
        varray_view(const varray_view<U>& other)
        {
            reset(other);
        }

        template <typename U, typename UAlloc, typename T_=T,
                  typename=detail::enable_if_const_t<T_>,
                  typename=detail::enable_if_convertible_t<U*,pointer>>
        varray_view(const varray<U, UAlloc>& other)
        {
            reset(other);
        }

        template <typename U, typename UAlloc,
                  typename=detail::enable_if_convertible_t<U*,pointer>>
        varray_view(varray<U, UAlloc>& other)
        {
            reset(other);
        }

        varray_view(const std::vector<idx_type>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset(len, ptr, layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        varray_view(const std::vector<U>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset(len, ptr, layout);
        }

        varray_view(const std::vector<idx_type>& len, pointer ptr, const std::vector<stride_type>& stride)
        {
            reset(len, ptr, stride);
        }

        template <typename U, typename V, typename=
                  detail::enable_if_t<std::is_integral<U>::value &&
                                      std::is_integral<V>::value>>
        varray_view(const std::vector<U>& len, pointer ptr, const std::vector<V>& stride)
        {
            reset(len, ptr, stride);
        }

        const varray_view& operator=(const varray_view& other) const
        {
            copy(other, *this);
            return *this;
        }

        template <typename U, typename=detail::enable_if_assignable_t<reference,U>>
        const varray_view& operator=(const varray_view<U>& other) const
        {
            copy(other, *this);
            return *this;
        }

        template <typename U, typename Alloc,
                  typename=detail::enable_if_assignable_t<reference,U>>
        const varray_view& operator=(const varray<U, Alloc>& other) const
        {
            copy(other.view(), *this);
            return *this;
        }

        const varray_view& operator=(const T& value) const
        {
            copy(value, *this);
            return *this;
        }

        void reset()
        {
            data_ = nullptr;
            len_.clear();
            stride_.clear();
        }

        void reset(varray_view&& other)
        {
            swap(other);
        }

        template <typename U>
        detail::enable_if_convertible_t<U*,pointer>
        reset(const varray_view<U>& other)
        {
            data_ = other.data();
            len_ = other.lengths();
            stride_ = other.strides();
        }

        template <typename U, typename UAlloc, typename T_=T,
                  typename=detail::enable_if_const_t<T_>>
        detail::enable_if_convertible_t<U*,pointer>
        reset(const varray<U, UAlloc>& other)
        {
            data_ = other.data();
            len_ = other.lengths();
            stride_ = other.strides();
        }

        template <typename U, typename UAlloc>
        detail::enable_if_convertible_t<U*,pointer>
        reset(varray<U, UAlloc>& other)
        {
            data_ = other.data();
            len_ = other.lengths();
            stride_ = other.strides();
        }

        void reset(const std::vector<idx_type>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset<idx_type>(len, ptr, layout);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::vector<U>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset(len, ptr, varray<typename std::remove_cv<T>::type>::default_strides(len, layout));
        }

        void reset(const std::vector<idx_type>& len, pointer ptr, const std::vector<stride_type>& stride)
        {
            reset<idx_type, stride_type>(len, ptr, stride);
        }

        template <typename U, typename V>
        detail::enable_if_t<std::is_integral<U>::value &&
                            std::is_integral<V>::value>
        reset(const std::vector<U>& len, pointer ptr, const std::vector<V>& stride)
        {
            MARRAY_ASSERT(len.size() > 0);
            MARRAY_ASSERT(len.size() == stride.size());
            data_ = ptr;
            len_.assign(len.begin(), len.end());
            stride_.assign(stride.begin(), stride.end());
        }

        void shift(unsigned dim, idx_type n)
        {
            MARRAY_ASSERT(dim < dimension());
            data_ += n*stride_[dim];
        }

        void shift_down(unsigned dim)
        {
            shift(dim, len_[dim]);
        }

        void shift_up(unsigned dim)
        {
            shift(dim, -len_[dim]);
        }

        varray_view<T> shifted(unsigned dim, idx_type n) const
        {
            MARRAY_ASSERT(dim < dimension());
            varray_view<T> r(*this);
            r.shift(dim, n);
            return r;
        }

        varray_view<T> shifted_down(unsigned dim) const
        {
            return shifted(dim, len_[dim]);
        }

        varray_view<T> shifted_up(unsigned dim) const
        {
            return shifted(dim, -len_[dim]);
        }

        void permute(const std::vector<unsigned>& perm)
        {
            permute<unsigned>(perm);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        permute(const std::vector<U>& perm)
        {
            MARRAY_ASSERT(perm.size() == dimension());

            std::vector<idx_type> len(len_);
            std::vector<stride_type> stride(stride_);

            for (unsigned i = 0;i < dimension();i++)
            {
                MARRAY_ASSERT(0 <= perm[i] && perm[i] < dimension());
                for (unsigned j = 0;j < i;j++) MARRAY_ASSERT(perm[i] != perm[j]);
            }

            for (unsigned i = 0;i < dimension();i++)
            {
                len_[i] = len[perm[i]];
                stride_[i] = stride[perm[i]];
            }
        }

        varray_view<T> permuted(const std::vector<unsigned>& perm) const
        {
            return permuted<unsigned>(perm);
        }

        template <typename U>
        detail::enable_if_integral_t<U,varray_view<T>>
        permuted(const std::vector<U>& perm) const
        {
            varray_view<T> r(*this);
            r.permute(perm);
            return r;
        }

        void lower(const std::vector<unsigned>& split)
        {
            lower<unsigned>(split);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        lower(const std::vector<U>& split)
        {
            MARRAY_ASSERT(split.size() < dimension());

            unsigned newdim = split.size()+1;
            for (unsigned i = 0;i < newdim-1;i++)
            {
                MARRAY_ASSERT(split[i] <= dimension());
                if (i != 0) MARRAY_ASSERT(split[i-1] <= split[i]);
            }

            std::vector<idx_type> len = len_;
            std::vector<stride_type> stride = stride_;

            for (unsigned i = 0;i < newdim;i++)
            {
                unsigned begin = (i == 0 ? 0 : split[i-1]);
                unsigned end = (i == newdim-1 ? dimension()-1 : split[i]-1);
                if (begin > end) continue;

                if (stride[begin] < stride[end])
                {
                    len_[i] = len[end];
                    stride_[i] = stride[begin];
                    for (auto j = begin;j < end;j++)
                    {
                        MARRAY_ASSERT(stride[j+1] == stride[j]*len[j]);
                        len_[i] *= len[j];
                    }
                }
                else
                {
                    len_[i] = len[end];
                    stride_[i] = stride[end];
                    for (auto j = begin;j < end;j++)
                    {
                        MARRAY_ASSERT(stride[j] == stride[j+1]*len[j+1]);
                        len_[i] *= len[j];
                    }
                }
            }

            len_.resize(newdim);
            stride_.resize(newdim);
        }

        varray_view<T> lowered(const std::vector<unsigned>& split) const
        {
            return lowered<unsigned>(split);
        }

        template <typename U>
        detail::enable_if_integral_t<U,varray_view<T>>
        lowered(const std::vector<U>& split) const
        {
            varray_view<T> r(*this);
            r.lower(split);
            return r;
        }

        void rotate_dim(unsigned dim, idx_type shift)
        {
            MArray::rotate_dim(*this, dim, shift);
        }

        void rotate(const std::vector<idx_type>& shift)
        {
            rotate<idx_type>(shift);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        rotate(const std::vector<U>& shift)
        {
            MARRAY_ASSERT(shift.size() == dimension());
            for (unsigned dim = 0;dim < dimension();dim++)
            {
                rotate_dim(dim, shift[dim]);
            }
        }

        const_reference cfront() const
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data_[0];
        }

        reference front() const
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data_[0];
        }

        varray_view<const T> cfront(unsigned dim) const
        {
            MARRAY_ASSERT(dimension() > 1);
            MARRAY_ASSERT(dim < dimension());
            MARRAY_ASSERT(len_[dim] > 0);

            std::vector<idx_type> len(dimension()-1);
            std::vector<stride_type> stride(dimension()-1);

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, dimension()-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, dimension()-dim-1, stride.begin()+dim);

            return {len, data_, stride};
        }

        varray_view<T> front(unsigned dim) const
        {
            MARRAY_ASSERT(dimension() > 1);
            MARRAY_ASSERT(dim < dimension());
            MARRAY_ASSERT(len_[dim] > 0);

            std::vector<idx_type> len(dimension()-1);
            std::vector<stride_type> stride(dimension()-1);

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, dimension()-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, dimension()-dim-1, stride.begin()+dim);

            return {len, data_, stride};
        }

        const_reference cback() const
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data_[(len_[0]-1)*stride_[0]];
        }

        reference back() const
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data_[(len_[0]-1)*stride_[0]];
        }

        varray_view<const T> cback(unsigned dim) const
        {
            MARRAY_ASSERT(dimension() > 1);
            MARRAY_ASSERT(dim < dimension());
            MARRAY_ASSERT(len_[dim] > 0);

            std::vector<idx_type> len(dimension()-1);
            std::vector<stride_type> stride(dimension()-1);

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, dimension()-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, dimension()-dim-1, stride.begin()+dim);

            return {len, data_+(len_[dim]-1)*stride_[dim], stride};
        }

        varray_view<T> back(unsigned dim) const
        {
            MARRAY_ASSERT(dimension() > 1);
            MARRAY_ASSERT(dim < dimension());
            MARRAY_ASSERT(len_[dim] > 0);

            std::vector<idx_type> len(dimension()-1);
            std::vector<stride_type> stride(dimension()-1);

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, dimension()-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, dimension()-dim-1, stride.begin()+dim);

            return {len, data_+(len_[dim]-1)*stride_[dim], stride};
        }

        template <typename... Args>
        detail::enable_if_t<detail::are_indices_or_slices<Args...>::value &&
                            !detail::are_convertible<idx_type, Args...>::value,
                            varray_view<T>>
        operator()(Args&&... args) const
        {
            MARRAY_ASSERT(sizeof...(Args) == dimension());

            pointer ptr = data();
            std::vector<idx_type> len;
            std::vector<stride_type> stride;

            get_slice<0>(ptr, len, stride, std::forward<Args>(args)...);

            return {len, ptr, stride};
        }

        template <typename... Args>
        detail::enable_if_t<detail::are_convertible<idx_type, Args...>::value,
                            reference>
        operator()(Args&&... args) const
        {
            MARRAY_ASSERT(sizeof...(Args) == dimension());
            pointer ptr = data();
            get_reference<0>(ptr, std::forward<Args>(args)...);
            return *ptr;
        }

        const_pointer cdata() const
        {
            return data_;
        }

        pointer data() const
        {
            return data_;
        }

        pointer data(pointer ptr)
        {
            std::swap(ptr, data_);
            return ptr;
        }

        idx_type length(unsigned dim) const
        {
            MARRAY_ASSERT(dim < dimension());
            return len_[dim];
        }

        idx_type length(unsigned dim, idx_type len)
        {
            MARRAY_ASSERT(dim < dimension());
            std::swap(len, len_[dim]);
            return len;
        }

        const std::vector<idx_type>& lengths() const
        {
            return len_;
        }

        stride_type stride(unsigned dim) const
        {
            MARRAY_ASSERT(dim < dimension());
            return stride_[dim];
        }

        stride_type stride(unsigned dim, stride_type stride)
        {
            MARRAY_ASSERT(dim < dimension());
            std::swap(stride, stride_[dim]);
            return stride;
        }

        const std::vector<stride_type>& strides() const
        {
            return stride_;
        }

        unsigned dimension() const
        {
            return static_cast<unsigned>(len_.size());
        }

        void swap(varray_view& other)
        {
            using std::swap;
            swap(data_, other.data_);
            swap(len_, other.len_);
            swap(stride_, other.stride_);
        }

        friend void swap(varray_view& a, varray_view& b)
        {
            a.swap(b);
        }
};

template <unsigned NDim, typename T>
marray_view<T, NDim> fix(const varray_view<T>& other)
{
    MARRAY_ASSERT(NDim == other.dimension());

    std::array<idx_type, NDim> len;
    std::array<stride_type, NDim> stride;
    std::copy_n(other.lengths().begin(), NDim, len.begin());
    std::copy_n(other.strides().begin(), NDim, stride.begin());

    return {len, other.data(), stride};
}

template <typename T, unsigned NDim>
varray_view<T> vary(const marray_view<T, NDim>& other)
{
    std::vector<idx_type> len{other.lengths().begin(), other.lengths().end()};
    std::vector<stride_type> stride{other.strides().begin(), other.strides().end()};
    return {len, other.data(), stride};
}

template <typename T, unsigned NDim, unsigned NIndexed, typename... Dims>
varray_view<T> vary(const marray_slice<T, NDim, NIndexed, Dims...>& other)
{
    return vary(other.view());
}

template <typename T, unsigned NDim, typename Alloc>
varray_view<const T> vary(const marray<T, NDim, Alloc>& other)
{
    std::vector<idx_type> len{other.lengths().begin(), other.lengths().end()};
    std::vector<stride_type> stride{other.strides().begin(), other.strides().end()};
    return {len, other.data(), stride};
}

template <typename T, unsigned NDim, typename Alloc>
varray_view<T> vary(marray<T, NDim, Alloc>& other)
{
    std::vector<idx_type> len{other.lengths().begin(), other.lengths().end()};
    std::vector<stride_type> stride{other.strides().begin(), other.strides().end()};
    return {len, other.data(), stride};
}

template <typename T>
varray_view<const T> cview(const varray_view<T>& x)
{
    return x;
}

template <typename T>
varray_view<T> view(const varray_view<T>& x)
{
    return x;
}

}

#endif
