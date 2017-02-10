#ifndef _MARRAY_MARRAY_VIEW_HPP_
#define _MARRAY_MARRAY_VIEW_HPP_

#include "marray.hpp"

namespace MArray
{

template <typename Type, unsigned NDim>
class marray_view
{
    static_assert(NDim > 0, "NDim must be positive");

    public:
        typedef Type value_type;
        typedef Type* pointer;
        typedef const Type* const_pointer;
        typedef Type& reference;
        typedef const Type& const_reference;

    protected:
        pointer data_ = nullptr;
        std::array<idx_type,NDim> len_ = {};
        std::array<stride_type,NDim> stride_ = {};

    public:
        marray_view() {}

        marray_view(const marray_view& other)
        {
            reset(other);
        }

        template <typename U, typename=detail::enable_if_convertible_t<U*,pointer>>
        marray_view(const marray_view<U, NDim>& other)
        {
            reset(other);
        }

        template <typename U, unsigned OldNDim, unsigned NIndexed, typename... Dims,
                  typename=detail::enable_if_convertible_t<U*,pointer>>
        marray_view(const marray_slice<U, OldNDim, NIndexed, Dims...>& other)
        {
            reset(other);
        }

        template <typename U, typename Alloc, typename T_=Type,
                  typename=detail::enable_if_const_t<T_>,
                  typename=detail::enable_if_convertible_t<const U*,pointer>>
        marray_view(const marray<U, NDim, Alloc>& other)
        {
            reset(other);
        }

        template <typename U, typename Alloc,
                  typename=detail::enable_if_convertible_t<U*,pointer>>
        marray_view(marray<U, NDim, Alloc>& other)
        {
            reset(other);
        }

        marray_view(const std::array<idx_type, NDim>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset(len, ptr, layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        marray_view(const std::array<U, NDim>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset(len, ptr, layout);
        }

        marray_view(const std::array<idx_type, NDim>& len, pointer ptr, const std::array<stride_type, NDim>& stride)
        {
            reset(len, ptr, stride);
        }

        template <typename U, typename V, typename=
                  detail::enable_if_t<std::is_integral<U>::value &&
                                      std::is_integral<V>::value>>
        marray_view(const std::array<U, NDim>& len, pointer ptr, const std::array<V, NDim>& stride)
        {
            reset(len, ptr, stride);
        }

        const marray_view& operator=(const marray_view& other) const
        {
            assign_expr(*this, other);
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        const marray_view& operator=(const Expression& other) const
        {
            assign_expr(*this, other);
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        const marray_view& operator+=(const Expression& other) const
        {
            *this = *this + other;
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        const marray_view& operator-=(const Expression& other) const
        {
            *this = *this - other;
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        const marray_view& operator*=(const Expression& other) const
        {
            *this = *this * other;
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        const marray_view& operator/=(const Expression& other) const
        {
            *this = *this / other;
            return *this;
        }

        void reset()
        {
            data_ = nullptr;
            len_.fill(0);
            stride_.fill(0);
        }

        template <typename U>
        detail::enable_if_convertible_t<U*,pointer>
        reset(const marray_view<U, NDim>& other)
        {
            data_ = other.data();
            len_ = other.lengths();
            stride_ = other.strides();
        }

        template <typename U, unsigned OldNDim, unsigned NIndexed, typename... Dims>
        detail::enable_if_convertible_t<U*,pointer>
        reset(const marray_slice<U, OldNDim, NIndexed, Dims...>& other)
        {
            reset(view(other));
        }

        template <typename U, typename Alloc, typename T_=Type,
                  typename=detail::enable_if_const_t<T_>>
        detail::enable_if_convertible_t<const U*,pointer>
        reset(const marray<U, NDim, Alloc>& other)
        {
            data_ = other.data();
            len_ = other.lengths();
            stride_ = other.strides();
        }

        template <typename U, typename Alloc>
        detail::enable_if_convertible_t<U*,pointer>
        reset(marray<U, NDim, Alloc>& other)
        {
            data_ = other.data();
            len_ = other.lengths();
            stride_ = other.strides();
        }

        void reset(const std::array<idx_type, NDim>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset<idx_type>(len, ptr, layout);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::array<U, NDim>& len, pointer ptr, layout layout = layout::DEFAULT)
        {
            reset(len, ptr, marray<typename std::remove_cv<Type>::type, NDim>::default_strides(len, layout));
        }

        void reset(const std::array<idx_type, NDim>& len, pointer ptr, const std::array<stride_type, NDim>& stride)
        {
            reset<idx_type, stride_type>(len, ptr, stride);
        }

        template <typename U, typename V>
        detail::enable_if_t<std::is_integral<U>::value &&
                            std::is_integral<V>::value>
        reset(const std::array<U, NDim>& len, pointer ptr, const std::array<V, NDim>& stride)
        {
            data_ = ptr;
            std::copy_n(len.begin(), NDim, len_.begin());
            std::copy_n(stride.begin(), NDim, stride_.begin());
        }

        template <unsigned Dim>
        void shift(idx_type n)
        {
            static_assert(Dim < NDim, "Dim out of range");
            data_ += n*stride_[Dim];
        }

        void shift(unsigned dim, idx_type n)
        {
            MARRAY_ASSERT(dim < NDim);
            data_ += n*stride_[dim];
        }

        template <unsigned Dim>
        void shift_down()
        {
            shift<Dim>(len_[Dim]);
        }

        void shift_down(unsigned dim)
        {
            shift(dim, len_[dim]);
        }

        template <unsigned Dim>
        void shift_up()
        {
            shift<Dim>(-len_[Dim]);
        }

        void shift_up(unsigned dim)
        {
            shift(dim, -len_[dim]);
        }

        template <unsigned Dim>
        marray_view<Type, NDim> shifted(idx_type n) const
        {
            marray_view<Type,NDim> r(*this);
            r.shift<Dim>(n);
            return r;
        }

        marray_view<Type, NDim> shifted(unsigned dim, idx_type n) const
        {
            marray_view<Type,NDim> r(*this);
            r.shift(dim, n);
            return r;
        }

        template <unsigned Dim>
        marray_view<Type,NDim> shifted_down() const
        {
            return shifted<Dim>(len_[Dim]);
        }

        marray_view<Type,NDim> shifted_down(unsigned dim) const
        {
            return shifted(dim, len_[dim]);
        }

        template <unsigned Dim>
        marray_view<Type,NDim> shifted_up() const
        {
            return shifted<Dim>(-len_[Dim]);
        }

        marray_view<Type,NDim> shifted_up(unsigned dim) const
        {
            return shifted(dim, -len_[dim]);
        }

        void permute(const std::array<unsigned, NDim>& perm)
        {
            permute<unsigned>(perm);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        permute(const std::array<U, NDim>& perm)
        {
            std::array<idx_type, NDim> len = len_;
            std::array<stride_type, NDim> stride = stride_;

            for (unsigned i = 0;i < NDim;i++)
            {
                MARRAY_ASSERT(0 <= perm[i] && perm[i] < NDim);
                for (unsigned j = 0;j < i;j++)
                    MARRAY_ASSERT(perm[i] != perm[j]);
            }

            for (unsigned i = 0;i < NDim;i++)
            {
                len_[i] = len[perm[i]];
                stride_[i] = stride[perm[i]];
            }
        }

        marray_view<Type,NDim> permuted(const std::array<unsigned, NDim>& perm) const
        {
            return permuted<unsigned>(perm);
        }

        template <typename U>
        detail::enable_if_integral_t<U,marray_view<Type,NDim>>
        permuted(const std::array<U, NDim>& perm) const
        {
            marray_view<Type,NDim> r(*this);
            r.permute(perm);
            return r;
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2>
        transpose()
        {
            permute({1, 0});
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2,marray_view<Type, NDim>>
        transposed() const
        {
            return permuted({1, 0});
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2,marray_view<Type, NDim>>
        T() const
        {
            return transposed();
        }

        template <size_t NewNDim>
        marray_view<Type, NewNDim> lowered(const std::array<unsigned, NewNDim-1>& split) const
        {
            return lowered<unsigned, NewNDim-1>(split);
        }

        template <typename U, size_t NSplit>
        detail::enable_if_integral_t<U,marray_view<Type, NSplit+1>>
        lowered(const std::array<U, NSplit>& split) const
        {
            MARRAY_ASSERT(NSplit < NDim);

            for (unsigned i = 0;i < NSplit;i++)
            {
                MARRAY_ASSERT(split[i] <= NDim);
                if (i != 0) MARRAY_ASSERT(split[i-1] <= split[i]);
            }

            std::array<idx_type, NSplit+1> newlen;
            std::array<stride_type, NSplit+1> newstride;

            for (unsigned i = 0;i <= NSplit;i++)
            {
                int begin = (i == 0 ? 0 : split[i-1]);
                int end = (i == NSplit ? NDim-1 : split[i]-1);
                if (begin > end) continue;

                if (stride_[begin] < stride_[end])
                {
                    newlen[i] = len_[end];
                    newstride[i] = stride_[begin];
                    for (int j = begin;j < end;j++)
                    {
                        MARRAY_ASSERT(stride_[j+1] == stride_[j]*len_[j]);
                        newlen[i] *= len_[j];
                    }
                }
                else
                {
                    newlen[i] = len_[end];
                    newstride[i] = stride_[end];
                    for (int j = begin;j < end;j++)
                    {
                        MARRAY_ASSERT(stride_[j] == stride_[j+1]*len_[j+1]);
                        newlen[i] *= len_[j];
                    }
                }
            }

            return {newlen, data_, newstride};
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data_[0];
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        front() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data_[0];
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return cfront();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        front() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return front();
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront(unsigned dim) const
        {
            MARRAY_ASSERT(dim == 0);
            return cfront();
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        front(unsigned dim) const
        {
            MARRAY_ASSERT(dim == 0);
            return front();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cfront() const
        {
            return cfront(Dim);
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front() const
        {
            return front(Dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cfront(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data_, stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data_, stride};
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data_[(len_[0]-1)*stride_[0]];
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, reference>
        back() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data_[(len_[0]-1)*stride_[0]];
        }

        template <unsigned Dim, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return cback();
        }

        template <unsigned Dim, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, reference>
        back() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return back();
        }

        template <unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback(unsigned dim) const
        {
            MARRAY_ASSERT(dim == 0);
            return cback();
        }

        template <unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, reference>
        back(unsigned dim) const
        {
            MARRAY_ASSERT(dim == 0);
            return back();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cback() const
        {
            return cback(Dim);
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back() const
        {
            return back(Dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cback(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data_ + (len_[dim]-1)*stride_[dim], stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data_ + (len_[dim]-1)*stride_[dim], stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        operator[](idx_type i) const
        {
            MARRAY_ASSERT(i < len_[0]);
            return data_[i*stride_[0]];
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_slice<Type, NDim, 1>>
        operator[](idx_type i) const
        {
            MARRAY_ASSERT(i < len_[0]);
            return {*this, i};
        }

        template <typename I>
        marray_slice<Type, NDim, 1, slice_dim>
        operator[](const range_t<I>& x) const
        {
            MARRAY_ASSERT(x.front() <= x.back());
            MARRAY_ASSERT(x.front() >= 0 && x.back() <= len_[0]);
            return {*this, x};
        }

        marray_slice<Type, NDim, 1, slice_dim>
        operator[](all_t) const
        {
            return {*this, range(len_[0])};
        }

        marray_slice<Type, NDim, 0, bcast_dim>
        operator[](bcast_t) const
        {
            return {*this};
        }

        template <typename Arg, typename=
            detail::enable_if_t<detail::is_index_or_slice<Arg>::value>>
        auto operator()(Arg&& arg) const ->
        decltype((*this)[std::forward<Arg>(arg)])
        {
            return (*this)[std::forward<Arg>(arg)];
        }

        template <typename Arg, typename... Args, typename=
            detail::enable_if_t<sizeof...(Args) &&
                detail::are_indices_or_slices<Arg, Args...>::value>>
        auto operator()(Arg&& arg, Args&&... args) const ->
        decltype((*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...))
        {
            return (*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...);
        }

        void rotate_dim(unsigned dim, idx_type shift)
        {
            MArray::rotate_dim(*this, dim, shift);
        }

        template <unsigned Dim>
        void rotate_dim(idx_type shift)
        {
            rotate_dim(Dim, shift);
        }

        void rotate(const std::array<idx_type, NDim>& shift)
        {
            rotate<idx_type>(shift);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        rotate(const std::array<U, NDim>& shift)
        {
            for (unsigned dim = 0;dim < NDim;dim++)
            {
                rotate_dim(dim, shift[dim]);
            }
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

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, idx_type>
        length() const
        {
            return len_[0];
        }

        template <unsigned Dim>
        idx_type length() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            return len_[Dim];
        }

        idx_type length(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            return len_[dim];
        }

        template <unsigned Dim>
        idx_type length(idx_type len)
        {
            static_assert(Dim < NDim, "Dim out of range");
            std::swap(len, len_[Dim]);
            return len;
        }

        idx_type length(unsigned dim, idx_type len)
        {
            MARRAY_ASSERT(dim < NDim);
            std::swap(len, len_[dim]);
            return len;
        }

        const std::array<idx_type, NDim>& lengths() const
        {
            return len_;
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, stride_type>
        stride() const
        {
            return stride_[0];
        }

        template <unsigned Dim>
        stride_type stride() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            return stride_[Dim];
        }

        stride_type stride(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            return stride_[dim];
        }

        template <unsigned Dim>
        stride_type stride(stride_type stride)
        {
            static_assert(Dim < NDim, "Dim out of range");
            std::swap(stride, stride_[Dim]);
            return stride;
        }

        stride_type stride(unsigned dim, stride_type stride)
        {
            MARRAY_ASSERT(dim < NDim);
            std::swap(stride, stride_[dim]);
            return stride;
        }

        const std::array<stride_type, NDim>& strides() const
        {
            return stride_;
        }

        static constexpr unsigned dimension()
        {
            return NDim;
        }

        void swap(marray_view& other)
        {
            using std::swap;
            swap(data_, other.data_);
            swap(len_, other.len_);
            swap(stride_, other.stride_);
        }

        friend void swap(marray_view& a, marray_view& b)
        {
            a.swap(b);
        }

        operator const marray_view<const Type, NDim>&() const
        {
            return reinterpret_cast<const marray_view<const Type, NDim>&>(*this);
        }

        operator marray_view<const Type, NDim>&()
        {
            return reinterpret_cast<marray_view<const Type, NDim>&>(*this);
        }
};

template <typename T, unsigned NDim>
marray_view<const T, NDim> cview(const marray_view<T, NDim>& x)
{
    return x;
}

template <typename T, unsigned NDim>
marray_view<T, NDim> view(const marray_view<T, NDim>& x)
{
    return x;
}

}

#endif
