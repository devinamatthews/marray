#ifndef _MARRAY_MARRAY_SLICE_HPP_
#define _MARRAY_MARRAY_SLICE_HPP_

#include "marray.hpp"

namespace MArray
{

/*
 * Represents a part of an array, where the first NIndexed-1 out of NDim
 * Dimensions have either been indexed into (i.e. a single value
 * specified for that index) or sliced (i.e. a range of values specified).
 * The parameter NSliced specifies how many indices were sliced. The
 * reference may be converted into an array view (of Dimension
 * NDim-NIndexed+1+NSliced) or further indexed, but may not be used to modify
 * data.
 */
template <typename Type, unsigned NDim, unsigned NIndexed, unsigned NSliced>
class marray_slice
{
    template <typename, unsigned, unsigned, unsigned> friend class marray_slice;

    public:
        typedef typename marray_view<Type, NDim>::value_type value_type;
        typedef typename marray_view<Type, NDim>::const_pointer const_pointer;
        typedef typename marray_view<Type, NDim>::pointer pointer;
        typedef typename marray_view<Type, NDim>::const_reference const_reference;
        typedef typename marray_view<Type, NDim>::reference reference;

    protected:
        static constexpr unsigned CurDim = NIndexed+NSliced-1;
        static constexpr unsigned NextDim = NIndexed+NSliced;
        static constexpr unsigned NewNDim = NDim-NIndexed;
        static constexpr unsigned DimsLeft = NDim-NIndexed-NSliced;
        static constexpr bool Const = std::is_const<Type>::value;

        const std::array<idx_type, NDim>& len_;
        const std::array<stride_type, NDim>& stride_;
        pointer data_;
        std::array<unsigned, NSliced> dims_;
        std::array<idx_type, NSliced> slice_len_;
        std::array<stride_type, NSliced> slice_stride_;

    public:
        marray_slice(const marray_slice& other) = default;

        template <typename Array,
                  unsigned N1=NIndexed, unsigned N2=NSliced,
                  typename=detail::enable_if_t<N1==1 && N2==0>>
        marray_slice(Array&& array, idx_type i)
        : len_(array.lengths()), stride_(array.strides()),
          data_(array.data() + i*stride_[CurDim]) {}

        template <typename Array, typename I,
                  unsigned N1=NIndexed, unsigned N2=NSliced,
                  typename=detail::enable_if_t<N1==0 && N2==1>>
        marray_slice(Array&& array, const range_t<I>& slice)
        : len_(array.lengths()), stride_(array.strides()),
          data_(array.data() + slice.front()*stride_[CurDim]),
          dims_{CurDim}, slice_len_{slice.size()},
          slice_stride_{slice.step()*stride_[CurDim]} {}

        marray_slice(const marray_slice<Type, NDim, NIndexed-1, NSliced>& parent, idx_type i)
        : len_(parent.len_), stride_(parent.stride_),
          data_(parent.data_ + i*parent.stride_[CurDim]),
          dims_(parent.dims_), slice_len_(parent.slice_len_),
          slice_stride_(parent.slice_stride_) {}

        template <typename I>
        marray_slice(const marray_slice<Type, NDim, NIndexed, NSliced-1>& parent,
                     const range_t<I>& slice)
        : len_(parent.len_), stride_(parent.stride_),
          data_(parent.data_ + slice.front()*parent.stride_[CurDim])
        {
            std::copy_n(parent.dims_.begin(), NSliced-1, dims_.begin());
            dims_.back() = CurDim;
            std::copy_n(parent.slice_len_.begin(), NSliced-1, slice_len_.begin());
            slice_len_.back() = slice.size();
            std::copy_n(parent.slice_stride_.begin(), NSliced-1, slice_stride_.begin());
            slice_stride_.back() = slice.step()*stride_[CurDim];
        }

        const marray_slice& operator=(const marray_slice& other) const
        {
            assign_expr(*this, other);
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        const marray_slice& operator=(const Expression& other) const
        {
            assign_expr(*this, other);
            return *this;
        }

        template <int N=NewNDim>
        detail::enable_if_t<N==1, reference>
        operator[](idx_type i) const
        {
            MARRAY_ASSERT(i >= 0 && i < len_[NextDim]);
            return data_[i*stride_[NextDim]];
        }

        template <int N=NewNDim>
        detail::enable_if_t<(N>1), marray_slice<Type, NDim, NIndexed+1, NSliced>>
        operator[](idx_type i) const
        {
            MARRAY_ASSERT(i >= 0 && i < len_[NextDim]);
            return {*this, i};
        }

        template <typename I>
        marray_slice<Type, NDim, NIndexed, NSliced+1>
        operator[](const range_t<I>& x) const
        {
            MARRAY_ASSERT(x.front() <= x.back());
            MARRAY_ASSERT(x.front() >= 0 && x.back() <= len_[NextDim]);
            return {*this, x};
        }

        marray_slice<Type, NDim, NIndexed, NSliced+1>
        operator[](all_t) const
        {
            return {*this, range(idx_type(), len_[NextDim])};
        }

        template <typename Arg, typename=
            detail::enable_if_t<DimsLeft == 1 &&
                                detail::is_index_or_slice<Arg>::value>>
        auto operator()(Arg&& arg) const ->
        decltype((*this)[std::forward<Arg>(arg)])
        {
            return (*this)[std::forward<Arg>(arg)];
        }

        template <typename Arg1, typename Arg2, typename... Args, typename=
            detail::enable_if_t<sizeof...(Args) == DimsLeft-2 &&
                                detail::are_indices_or_slices<Arg1, Arg2, Args...>::value>>
        auto operator()(Arg1&& arg1, Arg2&& arg2, Args&&... args) const ->
        decltype((*this)[std::forward<Arg1>(arg1)](std::forward<Arg2>(arg2), std::forward<Args>(args)...))
        {
            return (*this)[std::forward<Arg1>(arg1)](std::forward<Arg2>(arg2), std::forward<Args>(args)...);
        }

        const_pointer cdata() const
        {
            return data_;
        }

        pointer data() const
        {
            return data_;
        }

        idx_type slice_length(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NSliced);
            return slice_len_[dim];
        }

        template <unsigned Dim>
        idx_type slice_length() const
        {
            static_assert(Dim < NSliced, "Dim out of range");
            return slice_len_[Dim];
        }

        stride_type slice_stride(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NSliced);
            return slice_stride_[dim];
        }

        template <unsigned Dim>
        stride_type slice_stride() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            return slice_stride_[Dim];
        }

        idx_type base_length(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            return len_[dim];
        }

        template <unsigned Dim>
        idx_type base_length() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            return len_[Dim];
        }

        stride_type base_stride(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            return stride_[dim];
        }

        template <unsigned Dim>
        stride_type base_stride() const
        {
            static_assert(Dim < NSliced, "Dim out of range");
            return stride_[Dim];
        }

        marray_view<const Type, NewNDim> cview() const
        {
            std::array<idx_type, NewNDim> len;
            std::array<stride_type, NewNDim> stride;

            std::copy_n(slice_len_.begin(), NSliced, len.begin());
            std::copy_n(slice_stride_.begin(), NSliced, stride.begin());
            std::copy_n(len_.begin()+NextDim, DimsLeft, len.begin()+NSliced);
            std::copy_n(stride_.begin()+NextDim, DimsLeft, stride.begin()+NSliced);

            return {len, data_, stride};
        }

        marray_view<Type, NewNDim> view() const
        {
            std::array<idx_type, NewNDim> len;
            std::array<stride_type, NewNDim> stride;

            std::copy_n(slice_len_.begin(), NSliced, len.begin());
            std::copy_n(slice_stride_.begin(), NSliced, stride.begin());
            std::copy_n(len_.begin()+NextDim, DimsLeft, len.begin()+NSliced);
            std::copy_n(stride_.begin()+NextDim, DimsLeft, stride.begin()+NSliced);

            return {len, data_, stride};
        }

        friend marray_view<const Type, NewNDim> cview(const marray_slice& x)
        {
            return x.cview();
        }

        friend marray_view<Type, NewNDim> view(const marray_slice& x)
        {
            return x.view();
        }
};

}

#endif
