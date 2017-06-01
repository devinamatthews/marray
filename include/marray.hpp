#ifndef _MARRAY_MARRAY_HPP_
#define _MARRAY_MARRAY_HPP_

#include "utility.hpp"
#include "range.hpp"
#include "rotate.hpp"
#include "memory.hpp"

namespace MArray
{

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

template <typename Type, unsigned NDim, typename Allocator>
class marray : public marray_base<Type, NDim, true>
{
    static_assert(NDim > 0, "NDim must be positive");

    protected:
        typedef marray_base<Type, NDim, true> base;

    public:
        marray() {}

        marray(const marray& other)
        {
            reset(other);
        }

        marray(marray&& other)
        {
            reset(std::move(other));
        }

        template <typename U, bool Owner>
        marray(const marray_base<U, NDim, Owner>& other, layout layout = layout::DEFAULT)
        {
            reset(other, layout);
        }

        template <typename U, unsigned OldNDim, unsigned NIndexed, typename... Dims>
        marray(const marray_slice<U, OldNDim, NIndexed, Dims...>& other, layout layout = layout::DEFAULT)
        {
            reset(other, layout);
        }

        explicit marray(const std::array<idx_type, NDim>& len, const Type& val=Type(), layout layout = layout::DEFAULT)
        {
            reset(len, val, layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        explicit marray(const std::array<U, NDim>& len, const Type& val=Type(), layout layout = layout::DEFAULT)
        {
            reset(len, val, layout);
        }

        marray(const std::array<idx_type, NDim>& len, layout layout)
        {
            reset(len, Type(), layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        marray(const std::array<U, NDim>& len, layout layout)
        {
            reset(len, Type(), layout);
        }

        marray(const std::array<idx_type, NDim>& len, uninitialized_t, layout layout = layout::DEFAULT)
        {
            reset(len, uninitialized, layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        marray(const std::array<U, NDim>& len, uninitialized_t, layout layout = layout::DEFAULT)
        {
            reset(len, uninitialized, layout);
        }

        marray& operator=(const marray& other)
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

        void reset()
        {
            if (mem_)
            {
                mem_.destroy();
                mem_.deallocate();
            }

            len_.fill(0);
            stride_.fill(0);
            layout_ = layout::DEFAULT;
        }

        template <typename U>
        detail::enable_if_assignable_t<reference, U>
        reset(const marray_view<U, NDim>& other, layout layout = layout::DEFAULT)
        {
            if (std::is_scalar<Type>::value)
            {
                reset(other.lengths(), uninitialized, layout);
            }
            else
            {
                reset(other.lengths(), Type(), layout);
            }

            *this = other;
        }

        template <typename U, unsigned OldNDim, unsigned NIndexed, typename... Dims>
        detail::enable_if_assignable_t<reference, U>
        reset(const marray_slice<U, OldNDim, NIndexed, Dims...>& other, layout layout = layout::DEFAULT)
        {
            reset(other.view(), layout);
        }

        template <typename U, typename UAlloc>
        detail::enable_if_assignable_t<reference, U>
        reset(const marray<U, NDim, UAlloc>& other, layout layout = layout::DEFAULT)
        {
            reset(other.view(), layout);
        }

        void reset(marray&& other)
        {
            swap(other);
        }

        void reset(const std::array<idx_type, NDim>& len, const Type& val=Type(), layout layout = layout::DEFAULT)
        {
            reset<idx_type>(len, val, layout);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::array<U, NDim>& len, const Type& val=Type(), layout layout = layout::DEFAULT)
        {
            reset(len, uninitialized, layout);
            std::uninitialized_fill_n(mem_.data_, mem_.size_, val);
        }

        void reset(const std::array<idx_type, NDim>& len, layout layout)
        {
            reset<idx_type>(len, Type(), layout);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::array<U, NDim>& len, layout layout)
        {
            reset(len, Type(), layout);
        }

        void reset(const std::array<idx_type, NDim>& len, uninitialized_t, layout layout = layout::DEFAULT)
        {
            reset<idx_type>(len, uninitialized, layout);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::array<U, NDim>& len, uninitialized_t, layout layout = layout::DEFAULT)
        {
            reset();

            size_t size = std::accumulate(len.begin(), len.end(), size_type(1), std::multiplies<size_type>());
            mem_.allocate(size);

            layout_ = layout;
            std::copy_n(len.begin(), NDim, len_.begin());
            std::copy_n(default_strides(len, layout).begin(), NDim, stride_.begin());
        }

        void resize(const std::array<idx_type, NDim>& len, const Type& val=Type())
        {
            resize<idx_type>(len, val);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        resize(const std::array<U, NDim>& len, const Type& val=Type())
        {
            marray a(std::move(*this));
            reset(len, val, layout_);
            marray_view<Type, NDim> b(*this);

            /*
             * It is OK to change the geometry of 'a' even if it is not
             * a view since it is about to go out of scope.
             */
            for (unsigned i = 0;i < NDim;i++)
            {
                idx_type l = std::min(a.length(i), b.length(i));
                a.len_[i] = l;
                b.length(i, l);
            }

            b = a;
        }

        template <typename=void, unsigned N=NDim>
        typename std::enable_if<N==1>::type
        push_back(const Type& x)
        {
            resize({len_[0]+1});
            back() = x;
        }

        template <unsigned Dim, unsigned N=NDim>
        typename std::enable_if<N==1>::type
        push_back(const Type& x)
        {
            static_assert(Dim == 0, "Dim out of range");
            push_back(x);
        }

        template <unsigned N=NDim>
        typename std::enable_if<N==1>::type
        push_back(unsigned dim, const Type& x)
        {
            MARRAY_ASSERT(dim == 0);
            push_back(x);
        }

        template <unsigned Dim, typename U, unsigned N=NDim,
                  typename=detail::enable_if_assignable_t<reference, U>>
        void push_back(const marray_view<U, NDim-1>& x)
        {
            push_back(Dim, x);
        }

        template <typename U, unsigned N=NDim,
                  typename=detail::enable_if_assignable_t<reference, U>>
        void push_back(unsigned dim, const marray_view<U, NDim-1>& x)
        {
            MARRAY_ASSERT(dim < NDim);

            for (unsigned i = 0, j = 0;i < NDim;i++)
            {
                MARRAY_ASSERT(i == dim || len_[i] == x.length(j++));
            }

            std::array<idx_type, NDim> len = len_;
            len[dim]++;
            resize(len);
            back<NDim>(dim) = x;
        }

        template <typename=void, unsigned N=NDim>
        typename std::enable_if<N==1>::type
        pop_back()
        {
            resize({len_[0]-1});
        }

        template <unsigned Dim, unsigned N=NDim>
        void pop_back()
        {
            pop_back(Dim);
        }

        void pop_back(unsigned dim)
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim> len = len_;
            len[dim]--;
            resize(len);
        }

        marray_view<const Type, NDim> cview() const
        {
            return {len_, data(), stride_};
        }

        marray_view<const Type, NDim> view() const
        {
            return {len_, data(), stride_};
        }

        marray_view<Type, NDim> view()
        {
            return {len_, data(), stride_};
        }

        friend marray_view<const Type, NDim> cview(const marray& x)
        {
            return x.view();
        }

        friend marray_view<const Type, NDim> view(const marray& x)
        {
            return x.view();
        }

        friend marray_view<Type, NDim> view(marray& x)
        {
            return x.view();
        }

        marray_view<const Type, NDim> permuted(const std::array<unsigned, NDim>& perm) const
        {
            return permuted<unsigned>(perm);
        }

        marray_view<Type, NDim> permuted(const std::array<unsigned, NDim>& perm)
        {
            return permuted<unsigned>(perm);
        }

        template <typename U>
        detail::enable_if_integral_t<U, marray_view<const Type, NDim>>
        permuted(const std::array<U, NDim>& perm) const
        {
            marray_view<const Type, NDim> v = view();
            v.permute(perm);
            return v;
        }

        template <typename U>
        detail::enable_if_integral_t<U, marray_view<Type, NDim>>
        permuted(const std::array<U, NDim>& perm)
        {
            marray_view<Type, NDim> v = view();
            v.permute(perm);
            return v;
        }

        template <unsigned N=NDim>
        typename std::enable_if<N==2, marray_view<const Type, NDim>>::type
        transposed() const
        {
            marray_view<const Type, NDim> v = view();
            v.transpose();
            return v;
        }

        template <unsigned N=NDim>
        typename std::enable_if<N==2, marray_view<Type, NDim>>::type
        transposed()
        {
            marray_view<Type, NDim> v = view();
            v.transpose();
            return v;
        }

        template <unsigned N=NDim>
        typename std::enable_if<N==2, marray_view<const Type, NDim>>::type
        T() const
        {
            return transposed();
        }

        template <unsigned N=NDim>
        typename std::enable_if<N==2, marray_view<Type, NDim>>::type
        T()
        {
            return transposed();
        }

        template <size_t NewNDim>
        marray_view<const Type, NewNDim>
        lowered(const std::array<unsigned, NewNDim-1>& split) const
        {
            return view().lowered(split);
        }

        template <typename U, size_t NSplit>
        detail::enable_if_integral_t<U, marray_view<const Type, NSplit+1>>
        lowered(const std::array<U, NSplit>& split) const
        {
            return view().lowered(split);
        }

        template <size_t NewNDim>
        marray_view<Type, NewNDim>
        lowered(const std::array<unsigned, NewNDim-1>& split)
        {
            return view().lowered(split);
        }

        template <typename U, size_t NSplit>
        detail::enable_if_integral_t<U, marray_view<Type, NSplit+1>>
        lowered(const std::array<U, NSplit>& split)
        {
            return view().lowered(split);
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

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data()[0];
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        front() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data()[0];
        }

        template <typename=void, unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        front()
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data()[0];
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return cfront();
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        front() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return cfront();
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
            MARRAY_ASSERT(dim == 0);
            return cfront();
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        front(unsigned dim) const
        {
            MARRAY_ASSERT(dim == 0);
            return cfront();
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
            static_assert(Dim < NDim, "Dim out of range");
            MARRAY_ASSERT(len_[Dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), Dim, len.begin());
            std::copy_n(len_.begin()+Dim+1, NDim-Dim-1, len.begin()+Dim);
            std::copy_n(stride_.begin(), Dim, stride.begin());
            std::copy_n(stride_.begin()+Dim+1, NDim-Dim-1, stride.begin()+Dim);

            return {len, data(), stride};
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        front() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            MARRAY_ASSERT(len_[Dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), Dim, len.begin());
            std::copy_n(len_.begin()+Dim+1, NDim-Dim-1, len.begin()+Dim);
            std::copy_n(stride_.begin(), Dim, stride.begin());
            std::copy_n(stride_.begin()+Dim+1, NDim-Dim-1, stride.begin()+Dim);

            return {len, data(), stride};
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front()
        {
            static_assert(Dim < NDim, "Dim out of range");
            MARRAY_ASSERT(len_[Dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), Dim, len.begin());
            std::copy_n(len_.begin()+Dim+1, NDim-Dim-1, len.begin()+Dim);
            std::copy_n(stride_.begin(), Dim, stride.begin());
            std::copy_n(stride_.begin()+Dim+1, NDim-Dim-1, stride.begin()+Dim);

            return {len, data(), stride};
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

            return {len, data(), stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
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

            return {len, data(), stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front(unsigned dim)
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data(), stride};
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data()[(len_[0]-1)*stride_[0]];
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        back() const
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data()[(len_[0]-1)*stride_[0]];
        }

        template <typename=void, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, reference>
        back()
        {
            MARRAY_ASSERT(len_[0] > 0);
            return data()[(len_[0]-1)*stride_[0]];
        }

        template <unsigned Dim, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        cback() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return cback();
        }

        template <unsigned Dim, unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        back() const
        {
            static_assert(Dim == 0, "Dim out of range");
            return cback();
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
            MARRAY_ASSERT(dim == 0);
            return cback();
        }

        template <unsigned ndim_=NDim>
        detail::enable_if_t<ndim_==1, const_reference>
        back(unsigned dim) const
        {
            MARRAY_ASSERT(dim == 0);
            return cback();
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
            static_assert(Dim < NDim, "Dim out of range");
            MARRAY_ASSERT(len_[Dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), Dim, len.begin());
            std::copy_n(len_.begin()+Dim+1, NDim-Dim-1, len.begin()+Dim);
            std::copy_n(stride_.begin(), Dim, stride.begin());
            std::copy_n(stride_.begin()+Dim+1, NDim-Dim-1, stride.begin()+Dim);

            return {len, data() + (len_[Dim]-1)*stride_[Dim], stride};
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        back() const
        {
            static_assert(Dim < NDim, "Dim out of range");
            MARRAY_ASSERT(len_[Dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), Dim, len.begin());
            std::copy_n(len_.begin()+Dim+1, NDim-Dim-1, len.begin()+Dim);
            std::copy_n(stride_.begin(), Dim, stride.begin());
            std::copy_n(stride_.begin()+Dim+1, NDim-Dim-1, stride.begin()+Dim);

            return {len, data() + (len_[Dim]-1)*stride_[Dim], stride};
        }

        template <unsigned Dim, unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back()
        {
            static_assert(Dim < NDim, "Dim out of range");
            MARRAY_ASSERT(len_[Dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), Dim, len.begin());
            std::copy_n(len_.begin()+Dim+1, NDim-Dim-1, len.begin()+Dim);
            std::copy_n(stride_.begin(), Dim, stride.begin());
            std::copy_n(stride_.begin()+Dim+1, NDim-Dim-1, stride.begin()+Dim);

            return {len, data() + (len_[Dim]-1)*stride_[Dim], stride};
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

            return {len, data() + (len_[dim]-1)*stride_[dim], stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
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

            return {len, data() + (len_[dim]-1)*stride_[dim], stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back(unsigned dim)
        {
            MARRAY_ASSERT(dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<idx_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data() + (len_[dim]-1)*stride_[dim], stride};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, const_reference>
        operator[](idx_type i) const
        {
            MARRAY_ASSERT(i < len_[0]);
            return data()[i*stride_[0]];
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N==1, reference>
        operator[](idx_type i)
        {
            MARRAY_ASSERT(i < len_[0]);
            return data()[i*stride_[0]];
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_slice<const Type, NDim, 1>>
        operator[](idx_type i) const
        {
            MARRAY_ASSERT(i < len_[0]);
            return {*this, i};
        }

        template <unsigned N=NDim>
        detail::enable_if_t<N!=1, marray_slice<Type, NDim, 1>>
        operator[](idx_type i)
        {
            MARRAY_ASSERT(i < len_[0]);
            return {*this, i};
        }

        template <typename I>
        marray_slice<const Type, NDim, 1, slice_dim>
        operator[](const range_t<I>& x) const
        {
            MARRAY_ASSERT(x.front() <= x.back());
            MARRAY_ASSERT(x.front() >= 0 && x.back() <= len_[0]);
            return {*this, x};
        }

        template <typename I>
        marray_slice<Type, NDim, 1, slice_dim>
        operator[](const range_t<I>& x)
        {
            MARRAY_ASSERT(x.front() <= x.back());
            MARRAY_ASSERT(x.front() >= 0 && x.back() <= len_[0]);
            return {*this, x};
        }

        marray_slice<const Type, NDim, 1, slice_dim>
        operator[](all_t) const
        {
            return {*this, range(len_[0])};
        }

        marray_slice<Type, NDim, 1, slice_dim>
        operator[](all_t)
        {
            return {*this, range(len_[0])};
        }

        marray_slice<const Type, NDim, 0, bcast_dim>
        operator[](bcast_t) const
        {
            return {*this};
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

        const_pointer cdata() const
        {
            return mem_;
        }

        const_pointer data() const
        {
            return mem_;
        }

        pointer data()
        {
            return mem_;
        }

        template <typename=void, unsigned N=NDim>
        typename std::enable_if<N==1, idx_type>::type
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

        const std::array<idx_type, NDim>& lengths() const
        {
            return len_;
        }

        template <typename=void, unsigned N=NDim>
        typename std::enable_if<N==1, stride_type>::type
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

        const std::array<stride_type, NDim>& strides() const
        {
            return stride_;
        }

        static constexpr unsigned dimension()
        {
            return NDim;
        }

        void swap(marray& other)
        {
            using std::swap;
            swap(mem_, other.mem_);
            swap(layout_, other.layout_);
            base::swap(other);
        }

        friend void swap(marray& a, marray& b)
        {
            a.swap(b);
        }

    protected:
        memory<Type,Allocator> mem_;
        layout layout_ = layout::DEFAULT;
};

/*
 * Convenient names for 1- and 2-dimensional array types.
 */
template <typename Type> using row_view = marray_view<Type, 1>;
template <typename Type, typename Allocator=std::allocator<Type>> using row = marray<Type, 1, Allocator>;

template <typename Type> using matrix_view = marray_view<Type, 2>;
template <typename Type, typename Allocator=std::allocator<Type>> using matrix = marray<Type, 2, Allocator>;

}

#endif
