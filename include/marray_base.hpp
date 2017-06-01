#ifndef _MARRAY_MARRAY_BASE_HPP_
#define _MARRAY_MARRAY_BASE_HPP_

#include "utility.hpp"
#include "range.hpp"
#include "rotate.hpp"
#include "memory.hpp"

namespace MArray
{

template <typename Type, unsigned NDimLeft, typename Derived, bool Owner>
class marray_base;

template <typename Type, unsigned NDim, unsigned NIndexed, typename... Dims>
class marray_slice;

template <typename Type, unsigned NDim>
class marray_view;

template <typename Type, unsigned NDim, typename Allocator=std::allocator<Type>>
class marray;

template <typename Expr>
struct is_expression_arg_or_scalar;

}

#include "expression.hpp"
#include "marray_slice.hpp"
#include "marray_view.hpp"

namespace MArray
{

namespace detail
{

template <typename... Dims>
struct get_sliced_dims_helper;

template <>
struct get_sliced_dims_helper<>
{
    template <typename Lengths, typename Strides>
    get_sliced_dims_helper(Lengths lengths, Strides strides) {}
};

template <typename Dim, typename... Dims>
struct get_sliced_dims_helper<Dim, Dims...>
{
    template <typename Lengths, typename Strides>
    get_sliced_dims_helper(Lengths lengths, Strides strides,
                           const slice_dim& dim)
    {
        *lengths == dim.len;
        *strides == dim.stride;
    }

    template <typename Lengths, typename Strides>
    get_sliced_dims_helper(Lengths, Strides, const bcast_dim&)
    {
        static_assert(false, "Cannot take a view of a broadcasted array");
    }

    template <typename Lengths, typename Strides>
    get_sliced_dims_helper(Lengths lengths, Strides strides,
                           const Dim& dim, const Dims&... dims)
    {
        get_sliced_dim(lengths, strides, dim);
        get_sliced_dims_helper<Dims...>(++lengths, ++strides, dims...);
    }
};

template <typename Lengths, typename Strides,
          typename T, unsigned NDim, bool Owner, typename... Dims>
void get_sliced_dims(Lengths lengths, Strides strides,
                     const marray_base<T,NDim,Owner,Dims...>& array)
{
    get_sliced_dims_helper<Dims...>(lengths, strides, array.dims_);
}

}

template <typename Type, unsigned NDimLeft, typename Derived, bool Owner>
class marray_base
{
    protected:
        typedef typename std::conditional<Owner,const Type,Type>::type cType;
        typedef cType& cref;
        typedef cType* cptr;

        template <typename U>
        decltype(std::declval<U>().dims_) dims_helper_(U* u);
        std::tuple<> dims_helper(...);
        typedef decltype(dims_helper_((Derived*)0)) Dims;

        static constexpr unsigned NSliced = std::tuple_size<Dims>::value;
        static constexpr unsigned NDim = NDimLeft + NSliced;

    public:
        typedef Type value_type;
        typedef Type* pointer;
        typedef const Type* const_pointer;
        typedef Type& reference;
        typedef const Type& const_reference;

        template <unsigned N>
        static std::array<stride_type, N>
        default_strides(std::initializer_list<idx_type> len, layout layout = layout::DEFAULT)
        {
            return default_strides<N,decltype(len)>(len, layout);
        }

        template <unsigned N, typename U>
        static std::array<stride_type, N>
        default_strides(const U& len, layout layout = layout::DEFAULT)
        {
            //TODO: add alignment option

            std::array<stride_type, N> stride;

            if (N == 0) return stride;

            if (layout == layout::ROW_MAJOR)
            {
                stride[N-1] = 1;
                for (unsigned i = N;i --> 1;)
                {
                    stride[i-1] = stride[i]*len[i];
                }
            }
            else
            {
                stride[0] = 1;
                for (unsigned i = 1;i < N;i++)
                {
                    stride[i] = stride[i-1]*len[i-1];
                }
            }

            return stride;
        }

        marray_base& operator=(const marray_base& other)
        {
            assign_expr(*this, other);
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        marray& operator=(const Expression& other)
        {
            assign_expr(*this, other);
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        marray& operator+=(const Expression& other)
        {
            *this = *this + other;
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        marray& operator-=(const Expression& other)
        {
            *this = *this - other;
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        marray& operator*=(const Expression& other)
        {
            *this = *this * other;
            return *this;
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        marray& operator/=(const Expression& other)
        {
            *this = *this / other;
            return *this;
        }

        marray_view<const Type, NDim> cview() const
        {
            return const_cast<marray_base&>(*this).view();
        }

        marray_view<cType, NDim> view() const
        {
            return const_cast<marray_base&>(*this).view();
        }

        marray_view<Type, NDim> view()
        {
            std::array<idx_type,NDim> len;
            std::array<stride_type,NDim> stride;

            detail::get_sliced_dims(len.begin(), stride.begin(), *this);
            std::copy_n(lengths().begin(), NDimLeft, len.begin()+NSliced);
            std::copy_n(strides().begin(), NDimLeft, stride.begin()+NSliced);

            return {len, data(), stride};
        }

        friend marray_view<const Type, NDim> cview(const marray_base& x)
        {
            return x.view();
        }

        friend marray_view<cType, NDim> view(const marray_base& x)
        {
            return x.view();
        }

        friend marray_view<Type, NDim> view(marray_base& x)
        {
            return x.view();
        }

        marray_view<cType, NDim> shifted(std::initializer_list<idx_type> n) const
        {
            return const_cast<marray_base&>(*this).shifted(n);
        }

        marray_view<Type, NDim> shifted(std::initializer_list<idx_type> n)
        {
            marray_view<Type,NDim> r(*this);
            r.shift(n);
            return r;
        }

        template <unsigned Dim>
        marray_view<cType, NDim> shifted(idx_type n) const
        {
            return const_cast<marray_base&>(*this).shifted<Dim>(n);
        }

        template <unsigned Dim>
        marray_view<Type, NDim> shifted(idx_type n)
        {
            return shifted(Dim, n);
        }

        marray_view<cType, NDim> shifted(unsigned dim, idx_type n) const
        {
            return const_cast<marray_base&>(*this).shifted(dim, n);
        }

        marray_view<Type, NDim> shifted(unsigned dim, idx_type n)
        {
            marray_view<Type,NDim> r(*this);
            r.shift(dim, n);
            return r;
        }

        template <unsigned Dim>
        marray_view<cType,NDim> shifted_down() const
        {
            return const_cast<marray_base&>(*this).shifted_down<Dim>();
        }

        template <unsigned Dim>
        marray_view<Type,NDim> shifted_down()
        {
            return shifted_down(Dim);
        }

        marray_view<cType,NDim> shifted_down(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).shifted_down(dim);
        }

        marray_view<Type,NDim> shifted_down(unsigned dim)
        {
            return shifted(dim, length(dim));
        }

        template <unsigned Dim>
        marray_view<cType,NDim> shifted_up() const
        {
            return const_cast<marray_base&>(*this).shifted_up<Dim>();
        }

        template <unsigned Dim>
        marray_view<Type,NDim> shifted_up()
        {
            return shifted_up(Dim);
        }

        marray_view<cType,NDim> shifted_up(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).shifted_up(dim);
        }

        marray_view<Type,NDim> shifted_up(unsigned dim)
        {
            return shifted(dim, -length(dim));
        }

        marray_view<cType,NDim> permuted(std::initializer_list<unsigned> perm) const
        {
            return const_cast<marray_base&>(*this).permuted(perm);
        }

        marray_view<Type,NDim> permuted(std::initializer_list<unsigned> perm)
        {
            return permuted<decltype(perm)>(perm);
        }

        template <typename U>
        marray_view<cType,NDim> permuted(const U& perm) const
        {
            return const_cast<marray_base&>(*this).permuted<U>(perm);
        }

        template <typename U>
        marray_view<Type,NDim> permuted(const U& perm)
        {
            marray_view<Type,NDim> r(*this);
            r.permute(perm);
            return r;
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2,marray_view<cType, NDim>>
        transposed() const
        {
            return const_cast<marray_base&>(*this).transposed();
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2,marray_view<Type, NDim>>
        transposed()
        {
            return permuted({1, 0});
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2,marray_view<cType, NDim>>
        T() const
        {
            return const_cast<marray_base&>(*this).T();
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2,marray_view<Type, NDim>>
        T()
        {
            return transposed();
        }

        template <size_t NewNDim>
        marray_view<cType, NewNDim> lowered(std::initializer_list<unsigned> split) const
        {
            return const_cast<marray_base&>(*this).lowered<NewNDim>(split);
        }

        template <size_t NewNDim>
        marray_view<Type, NewNDim> lowered(std::initializer_list<unsigned> split)
        {
            return lowered<NewNDim,decltype(split)>(split);
        }

        template <size_t NewNDim, typename U>
        marray_view<cType, NewNDim> lowered(const U& split) const
        {
            return const_cast<marray_base&>(*this).lowered<NewNDim,U>(split);
        }

        template <size_t NewNDim, typename U>
        marray_view<Type, NewNDim> lowered(const U& split)
        {
            constexpr size_t NSplit = NewNDim-1;
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

                if (stride(begin) < stride(end))
                {
                    newlen[i] = length(end);
                    newstride[i] = stride(begin);
                    for (int j = begin;j < end;j++)
                    {
                        MARRAY_ASSERT(stride(j+1) == stride(j)*length(j));
                        newlen[i] *= length(j);
                    }
                }
                else
                {
                    newlen[i] = length(end);
                    newstride[i] = stride(end);
                    for (int j = begin;j < end;j++)
                    {
                        MARRAY_ASSERT(stride(j) == stride(j+1)*length(j+1));
                        newlen[i] *= length(j);
                    }
                }
            }

            return {newlen, data(), newstride};
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront() const
        {
            return const_cast<marray_base&>(*this).front();
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, cref>
        front() const
        {
            return const_cast<marray_base&>(*this).front();
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        front()
        {
            MARRAY_ASSERT(length(0) > 0);
            return data()[0];
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N==1, cref>
        front() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        front()
        {
            static_assert(Dim == 0, "Dim out of range");
            return front();
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, cref>
        front(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        front(unsigned dim)
        {
            MARRAY_ASSERT(dim == 0);
            return front();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cfront() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<cType, NDim-1>>
        front() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front()
        {
            return front<N>(Dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cfront(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<cType, NDim-1>>
        front(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front(unsigned dim)
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(lengths().begin(), dim, len.begin());
            std::copy_n(lengths().begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(strides().begin(), dim, stride.begin());
            std::copy_n(strides().begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data(), stride};
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback() const
        {
            return const_cast<marray_base&>(*this).back();
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, cref>
        back() const
        {
            return const_cast<marray_base&>(*this).back();
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, reference>
        back()
        {
            MARRAY_ASSERT(length(0) > 0);
            return data()[(length(0)-1)*stride(0)];
        }

        template <unsigned Dim, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <unsigned Dim, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, cref>
        back() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <unsigned Dim, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, reference>
        back()
        {
            static_assert(Dim == 0, "Dim out of range");
            return back();
        }

        template <unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, cref>
        back(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, reference>
        back(unsigned dim)
        {
            MARRAY_ASSERT(dim == 0);
            return back();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cback() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<cType, NDim-1>>
        back() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back()
        {
            return back<N>(Dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cback(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<cType, NDim-1>>
        back(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back(unsigned dim)
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(length(dim) > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(lengths().begin(), dim, len.begin());
            std::copy_n(lengths().begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(strides().begin(), dim, stride.begin());
            std::copy_n(strides().begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data() + (length(dim)-1)*stride(dim), stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, cref>
        operator[](idx_type i) const
        {
            return const_cast<marray_base&>(*this)[i];
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        operator[](idx_type i)
        {
            MARRAY_ASSERT(i < length(0));
            return data()[i*stride(0)];
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_slice<cType, NDim, 1>>
        operator[](idx_type i) const
        {
            return const_cast<marray_base&>(*this)[i];
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_slice<Type, NDim, 1>>
        operator[](idx_type i)
        {
            MARRAY_ASSERT(i < length(0));
            return {*this, i};
        }

        template <typename I>
        marray_slice<cType, NDim, 1, slice_dim>
        operator[](const range_t<I>& x) const
        {
            return const_cast<marray_base&>(*this)[x];
        }

        template <typename I>
        marray_slice<Type, NDim, 1, slice_dim>
        operator[](const range_t<I>& x)
        {
            MARRAY_ASSERT(x.front() <= x.back());
            MARRAY_ASSERT(x.front() >= 0 && x.back() <= length(0));
            return {*this, x};
        }

        marray_slice<cType, NDim, 1, slice_dim>
        operator[](all_t) const
        {
            return const_cast<marray_base&>(*this)[slice::all];
        }

        marray_slice<Type, NDim, 1, slice_dim>
        operator[](all_t)
        {
            return {*this, range(length(0))};
        }

        marray_slice<cType, NDim, 0, bcast_dim>
        operator[](bcast_t) const
        {
            return const_cast<marray_base&>(*this)[slice::bcast];
        }

        marray_slice<Type, NDim, 0, bcast_dim>
        operator[](bcast_t)
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

        template <typename Arg, typename=
            detail::enable_if_t<detail::is_index_or_slice<Arg>::value>>
        auto operator()(Arg&& arg) ->
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

        template <typename Arg, typename... Args, typename=
            detail::enable_if_t<sizeof...(Args) &&
                detail::are_indices_or_slices<Arg, Args...>::value>>
        auto operator()(Arg&& arg, Args&&... args) ->
        decltype((*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...))
        {
            return (*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...);
        }

        template <bool O=Owner>
        typename std::enable_if<O>::type
        rotate_dim(unsigned dim, idx_type shift)
        {
            rotate_dim<false>(dim, shift);
        }

        template <bool O=Owner>
        typename std::enable_if<!O>::type
        rotate_dim(unsigned dim, idx_type shift) const
        {
            MArray::rotate_dim(*this, dim, shift);
        }

        template <unsigned Dim, bool O=Owner>
        typename std::enable_if<O>::type
        rotate_dim(idx_type shift)
        {
            rotate_dim<Dim,false>(shift);
        }

        template <unsigned Dim, bool O=Owner>
        typename std::enable_if<!O>::type
        rotate_dim(idx_type shift) const
        {
            rotate_dim(Dim, shift);
        }

        template <bool O=Owner>
        typename std::enable_if<O>::type
        rotate(std::initializer_list<idx_type> shift)
        {
            rotate<false>(shift);
        }

        template <bool O=Owner>
        typename std::enable_if<!O>::type
        rotate(std::initializer_list<idx_type> shift) const
        {
            rotate<decltype(shift)>(shift);
        }

        template <typename U, bool O=Owner>
        typename std::enable_if<O>::type
        rotate(const U& shift)
        {
            rotate<U,false>(shift);
        }

        template <typename U, bool O=Owner>
        typename std::enable_if<!O>::type
        rotate(const U& shift) const
        {
            for (unsigned dim = 0;dim < NDim;dim++)
            {
                rotate_dim(dim, shift[dim]);
            }
        }

        const_pointer cdata() const
        {
            return const_cast<marray_base&>(*this).data();
        }

        cptr data() const
        {
            return const_cast<marray_base&>(*this).data();
        }

        pointer data()
        {
            return derived().data_;
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, idx_type>
        length() const
        {
            return length(0);
        }

        template <unsigned Dim>
        idx_type length() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            return length(Dim);
        }

        idx_type length(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            return lengths()[dim];
        }

        const std::array<idx_type, NDim>& lengths() const
        {
            return derived().len_;
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, stride_type>
        stride() const
        {
            return stride(0);
        }

        template <unsigned Dim>
        stride_type stride() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            return stride(Dim);
        }

        stride_type stride(unsigned dim) const
        {
            MARRAY_ASSERT(dim < NDim);
            return strides()[dim];
        }

        const std::array<stride_type, NDim>& strides() const
        {
            return derived().stride_;
        }

        static constexpr unsigned dimension()
        {
            return NDim;
        }

    protected:
        void shift(std::initializer_list<idx_type> n)
        {
            MARRAY_ASSERT(n.size() == NDim);
            for (unsigned dim = 0;dim < NDim;dim++)
            {
                shift(dim, n.begin()[dim]);
            }
        }

        template <unsigned Dim>
        void shift(idx_type n)
        {
            static_assert(Dim < NDim, "Dim out of range");
            shift(Dim, n);
        }

        void shift(unsigned dim, idx_type n)
        {
            MARRAY_ASSERT(dim < NDim);
            derived().data_ += n*derived().stride_[dim];
        }

        template <unsigned Dim>
        void shift_down()
        {
            shift_down(Dim);
        }

        void shift_down(unsigned dim)
        {
            shift(dim, derived().len_[dim]);
        }

        template <unsigned Dim>
        void shift_up()
        {
            shift_up(Dim);
        }

        void shift_up(unsigned dim)
        {
            shift(dim, -derived().len_[dim]);
        }

        void permute(const std::array<unsigned, NDim>& perm)
        {
            permute<decltype(perm)>(perm);
        }

        template <typename U>
        void permute(const U& perm)
        {
            std::array<idx_type, NDim> len = derived().len_;
            std::array<stride_type, NDim> stride = derived().stride_;

            for (unsigned i = 0;i < NDim;i++)
            {
                MARRAY_ASSERT(0 <= perm[i] && perm[i] < NDim);
                for (unsigned j = 0;j < i;j++)
                    MARRAY_ASSERT(perm[i] != perm[j]);
            }

            for (unsigned i = 0;i < NDim;i++)
            {
                derived().len_[i] = len[perm[i]];
                derived().stride_[i] = stride[perm[i]];
            }
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==2> transpose()
        {
            permute({1, 0});
        }

        pointer data(pointer ptr)
        {
            std::swap(ptr, derived().data_);
            return ptr;
        }

        template <unsigned Dim>
        idx_type length(idx_type len)
        {
            static_assert(Dim < NDim, "Dim out of range");
            return length(Dim, len);
        }

        idx_type length(unsigned dim, idx_type len)
        {
            MARRAY_ASSERT(dim < NDim);
            std::swap(len, derived().len_[dim]);
            return len;
        }

        template <unsigned Dim>
        stride_type stride(stride_type s)
        {
            static_assert(Dim < NDim, "Dim out of range");
            return stride(Dim, s);
        }

        stride_type stride(unsigned dim, stride_type s)
        {
            MARRAY_ASSERT(dim < NDim);
            std::swap(s, derived().stride_[dim]);
            return s;
        }

        const Derived& derived() const
        {
            return static_cast<const Derived&>(*this);
        }

        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }
};

template <typename Type, unsigned NDim>
class marray_data
{
    protected:
        Type* data_;
        std::array<idx_type,NDim> len_;
        std::array<idx_type,NDim> stride_;
};

template <typename Type, unsigned NDim>
class marray_ref : public marray_base<Type, NDim, marray_ref, false>
{
    protected:
        const marray_data<Type, NDim>& data_;
};

}

#endif
