#ifndef _MARRAY_VARRAY_HPP_
#define _MARRAY_VARRAY_HPP_

#include "marray.hpp"

namespace MArray
{

template <typename T>
class varray_view;

template <typename T, typename Allocator=std::allocator<T>>
class varray;

}

#include "varray_view.hpp"

namespace MArray
{

template <typename T, typename U>
detail::enable_if_assignable_t<U&,T>
copy(const varray_view<T>& a, const varray_view<U>& b);

template <typename T, typename U>
detail::enable_if_assignable_t<U&,T>
copy(const T& a, const varray_view<U>& b);

template <typename T, typename Allocator>
class varray
{
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;

    protected:
        struct alloc_s_ : Allocator
        {
            typedef std::allocator_traits<Allocator> traits_;

            pointer data_ = nullptr;
            size_t size_ = 0;

            alloc_s_() {}

            alloc_s_(const Allocator& alloc) : Allocator(alloc) {}

            void allocate(size_t size)
            {
                data_ = traits_::allocate(*this, size);
                size_ = size;
            }

            void deallocate()
            {
                traits_::deallocate(*this, data_, size_);
                data_ = nullptr;
                size_ = 0;
            }

            void destroy()
            {
                for (size_t i = 0;i < size_;i++)
                {
                    traits_::destroy(*this, data_+i);
                }
            }

            operator const_pointer() const { return data_; }

            operator pointer() { return data_; }

            explicit operator bool() const { return data_; }
        } alloc_;
        std::vector<idx_type> len_;
        std::vector<stride_type> stride_;
        layout layout_ = layout::DEFAULT;

        template <unsigned Dim, typename Ptr>
        void get_slice(Ptr& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride) const {}

        template <unsigned Dim, typename Ptr, typename... Args>
        void get_slice(Ptr& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride,
                       idx_type arg, Args&&... args) const
        {
            MARRAY_ASSERT(arg >= 0 && arg < len_[Dim]);
            ptr += arg*stride_[Dim];
            get_slice<Dim+1>(ptr, len, stride, std::forward<Args>(args)...);
        }

        template <unsigned Dim, typename Ptr, typename I, typename... Args>
        void get_slice(Ptr& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride,
                       const range_t<I>& arg, Args&&... args) const
        {
            MARRAY_ASSERT(arg.front() <= arg.back());
            MARRAY_ASSERT(arg.front() >= 0 && arg.back() <= len_[Dim]);
            ptr += arg.front()*stride_[Dim];
            len.push_back(arg.size());
            stride.push_back(arg.step()*stride_[Dim]);
            get_slice<Dim+1>(ptr, len, stride, std::forward<Args>(args)...);
        }

        template <unsigned Dim, typename Ptr, typename... Args>
        void get_slice(Ptr& ptr, std::vector<idx_type>& len,
                       std::vector<stride_type>& stride,
                       all_t, Args&&... args) const
        {
            len.push_back(len_[Dim]);
            stride.push_back(stride_[Dim]);
            get_slice<Dim+1>(ptr, len, stride, std::forward<Args>(args)...);
        }

        template <unsigned Dim, typename Ptr>
        void get_reference(Ptr& ptr) const {}

        template <unsigned Dim, typename Ptr, typename... Args>
        void get_reference(Ptr& ptr, idx_type arg, Args&&... args) const
        {
            MARRAY_ASSERT(arg >= 0 && arg < len_[Dim]);
            ptr += arg*stride_[Dim];
            get_reference<Dim+1>(ptr, std::forward<Args>(args)...);
        }

    public:
        static std::vector<stride_type> default_strides(const std::vector<idx_type>& len, layout layout = layout::DEFAULT)
        {
            return default_strides<idx_type>(len, layout);
        }

        template <typename U>
        static detail::enable_if_integral_t<U,std::vector<stride_type>>
        default_strides(const std::vector<U>& len, layout layout = layout::DEFAULT)
        {
            std::vector<stride_type> stride(len.size());

            if (stride.empty()) return stride;

            auto ndim = len.size();
            if (layout == layout::ROW_MAJOR)
            {
                stride[ndim-1] = 1;
                for (auto i = ndim;i --> 1;)
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

        varray() {}

        varray(const varray& other)
        {
            reset(other);
        }

        varray(varray&& other)
        {
            reset(std::move(other));
        }

        template <typename U, typename=detail::enable_if_assignable_t<reference,U>>
        varray(const varray_view<U>& other, layout layout = layout::DEFAULT)
        {
            reset(other, layout);
        }

        template <typename U, typename UAlloc,
                  typename=detail::enable_if_assignable_t<reference,U>>
        varray(const varray<U, UAlloc>& other, layout layout = layout::DEFAULT)
        {
            reset(other, layout);
        }

        explicit varray(const std::vector<idx_type>& len, const T& val=T(), layout layout = layout::DEFAULT)
        {
            reset(len, val, layout);
        }

        varray(const std::vector<idx_type>& len, layout layout)
        {
            reset(len, T(), layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        explicit varray(const std::vector<U>& len, const T& val=T(), layout layout = layout::DEFAULT)
        {
            reset(len, val, layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        varray(const std::vector<U>& len, layout layout)
        {
            reset(len, T(), layout);
        }

        varray(const std::vector<idx_type>& len, uninitialized_t u, layout layout = layout::DEFAULT)
        {
            reset(len, u, layout);
        }

        template <typename U, typename=detail::enable_if_integral_t<U>>
        varray(const std::vector<U>& len, uninitialized_t u, layout layout = layout::DEFAULT)
        {
            reset(len, u, layout);
        }

        ~varray()
        {
            reset();
        }

        const varray& operator=(const varray& other)
        {
            copy(other.view(), view());
            return *this;
        }

        template <typename U, typename=detail::enable_if_assignable_t<reference,U>>
        const varray& operator=(const varray_view<U>& other)
        {
            copy(other, view());
            return *this;
        }

        template <typename U, typename UAlloc,
                  typename=detail::enable_if_assignable_t<reference,U>>
        const varray& operator=(const varray<U, UAlloc>& other)
        {
            copy(other.view(), view());
            return *this;
        }

        template <typename U, typename=detail::enable_if_assignable_t<reference,U>>
        const varray& operator=(const U& value)
        {
            copy(value, view());
            return *this;
        }

        void reset()
        {
            if (alloc_)
            {
                alloc_.destroy();
                alloc_.deallocate();
            }

            len_.clear();
            stride_.clear();
            layout_ = layout::DEFAULT;
        }

        void reset(varray&& other)
        {
            swap(other);
        }

        template <typename U>
        detail::enable_if_assignable_t<reference,U>
        reset(const varray_view<U>& other, layout layout = layout::DEFAULT)
        {
            if (std::is_scalar<T>::value)
            {
                reset(other.lengths(), uninitialized, layout);
            }
            else
            {
                reset(other.lengths(), T(), layout);
            }

            *this = other;
        }

        template <typename U, typename UAlloc>
        detail::enable_if_assignable_t<reference,U>
        reset(const varray<U, UAlloc>& other, layout layout = layout::DEFAULT)
        {
            if (std::is_scalar<T>::value)
            {
                reset(other.lengths(), uninitialized, layout);
            }
            else
            {
                reset(other.lengths(), T(), layout);
            }

            *this = other;
        }

        void reset(const std::vector<idx_type>& len, const T& val=T(), layout layout = layout::DEFAULT)
        {
            reset<idx_type>(len, val, layout);
        }

        void reset(const std::vector<idx_type>& len, layout layout)
        {
            reset<idx_type>(len, T(), layout);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::vector<U>& len, const T& val=T(), layout layout = layout::DEFAULT)
        {
            reset(len, uninitialized, layout);
            std::uninitialized_fill_n(alloc_.data_, alloc_.size_, val);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::vector<U>& len, layout layout)
        {
            reset(len, T(), layout);
        }

        void reset(const std::vector<idx_type>& len, uninitialized_t, layout layout = layout::DEFAULT)
        {
            reset<idx_type>(len, uninitialized, layout);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        reset(const std::vector<U>& len, uninitialized_t, layout layout = layout::DEFAULT)
        {
            MARRAY_ASSERT(len.size() > 0);

            reset();

            size_t size = std::accumulate(len.begin(), len.end(), size_type(1), std::multiplies<size_type>());
            alloc_.allocate(size);

            layout_ = layout;
            len_.assign(len.begin(), len.end());
            stride_ = default_strides(len, layout);
        }

        void resize(const std::vector<idx_type>& len, const T& val=T())
        {
            resize<idx_type>(len, val);
        }

        template <typename U>
        detail::enable_if_integral_t<U>
        resize(const std::vector<U>& len, const T& val=T())
        {
            MARRAY_ASSERT(len.size() == dimension());

            varray a(std::move(*this));
            reset(len, val, layout_);
            auto b = view();

            /*
             * It is OK to change the geometry of 'a' even if it is not
             * a view since it is about to go out of scope.
             */
            for (unsigned i = 0;i < dimension();i++)
            {
                idx_type len = std::min(a.length(i), b.length(i));
                a.len_[i] = len;
                b.length(i, len);
            }

            b = a;
        }

        void push_back(const T& x)
        {
            MARRAY_ASSERT(dimension() == 1);
            resize({len_[0]+1});
            back() = x;
        }

        template <typename U>
        detail::enable_if_assignable_t<reference,U>
        push_back(unsigned dim, const varray_view<U>& x)
        {
            MARRAY_ASSERT(x.dimension()+1 == dimension());
            MARRAY_ASSERT(dim < dimension());

            for (unsigned i = 0, j = 0;i < dimension();i++)
            {
                MARRAY_ASSERT(i == dim || len_[i] == x.length(j++));
            }

            std::vector<idx_type> len = len_;
            len[dim]++;
            resize(len);
            back(dim) = x;
        }

        void pop_back()
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            resize({len_[0]-1});
        }

        void pop_back(unsigned dim)
        {
            MARRAY_ASSERT(dim < dimension());
            MARRAY_ASSERT(len_[dim] > 0);

            std::vector<idx_type> len = len_;
            len[dim]--;
            resize(len);
        }

        varray_view<const T> cview() const
        {
            return {len_, data(), stride_};
        }

        varray_view<const T> view() const
        {
            return {len_, data(), stride_};
        }

        varray_view<T> view()
        {
            return {len_, data(), stride_};
        }

        friend varray_view<const T> cview(const varray& x)
        {
            return x.view();
        }

        friend varray_view<const T> view(const varray& x)
        {
            return x.view();
        }

        friend varray_view<T> view(varray& x)
        {
            return x.view();
        }

        varray_view<const T> permuted(const std::vector<unsigned>& perm) const
        {
            return permuted<unsigned>(perm);
        }

        varray_view<T> permuted(const std::vector<unsigned>& perm)
        {
            return permuted<unsigned>(perm);
        }

        template <typename U>
        detail::enable_if_integral_t<U,varray_view<const T>>
        permuted(const std::vector<U>& perm) const
        {
            varray_view<const T> v = view();
            v.permute(perm);
            return v;
        }

        template <typename U>
        detail::enable_if_integral_t<U,varray_view<T>>
        permuted(const std::vector<U>& perm)
        {
            varray_view<T> v = view();
            v.permute(perm);
            return v;
        }

        varray_view<const T> lowered(const std::vector<unsigned>& split) const
        {
            return lowered<unsigned>(split);
        }

        varray_view<T> lowered(const std::vector<unsigned>& split)
        {
            return lowered<unsigned>(split);
        }

        template <typename U>
        detail::enable_if_integral_t<U,varray_view<const T>>
        lowered(const std::vector<U>& split) const
        {
            varray_view<const T> v = view();
            v.lower(split);
            return v;
        }

        template <typename U>
        detail::enable_if_integral_t<U,varray_view<T>>
        lowered(const std::vector<U>& split)
        {
            varray_view<T> v = view();
            v.lower(split);
            return v;
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
            return data()[0];
        }

        const_reference front() const
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data()[0];
        }

        reference front()
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data()[0];
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

            return {len, data(), stride};
        }

        varray_view<const T> front(unsigned dim) const
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

            return {len, data(), stride};
        }

        varray_view<T> front(unsigned dim)
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

            return {len, data(), stride};
        }

        const_reference cback() const
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data()[(len_[0]-1)*stride_[0]];
        }

        const_reference back() const
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data()[(len_[0]-1)*stride_[0]];
        }

        reference back()
        {
            MARRAY_ASSERT(dimension() == 1);
            MARRAY_ASSERT(len_[0] > 0);
            return data()[(len_[0]-1)*stride_[0]];
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

            return {len, data()+(len_[dim]-1)*stride_[dim], stride};
        }

        varray_view<const T> back(unsigned dim) const
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

            return {len, data()+(len_[dim]-1)*stride_[dim], stride};
        }

        varray_view<T> back(unsigned dim)
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

            return {len, data()+(len_[dim]-1)*stride_[dim], stride};
        }

        template <typename... Args>
        detail::enable_if_t<detail::are_indices_or_slices<Args...>::value &&
                            !detail::are_convertible<idx_type, Args...>::value,
                            varray_view<const T>>
        operator()(Args&&... args) const
        {
            MARRAY_ASSERT(sizeof...(Args) == dimension());

            const_pointer ptr = data();
            std::vector<idx_type> len;
            std::vector<stride_type> stride;

            get_slice<0>(ptr, len, stride, std::forward<Args>(args)...);

            return {len, ptr, stride};
        }

        template <typename... Args>
        detail::enable_if_t<detail::are_indices_or_slices<Args...>::value &&
                            !detail::are_convertible<idx_type, Args...>::value,
                            varray_view<T>>
        operator()(Args&&... args)
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
                            const_reference>
        operator()(Args&&... args) const
        {
            MARRAY_ASSERT(sizeof...(Args) == dimension());
            const_pointer ptr = data();
            get_reference<0>(ptr, std::forward<Args>(args)...);
            return *ptr;
        }

        template <typename... Args>
        detail::enable_if_t<detail::are_convertible<idx_type, Args...>::value,
                            reference>
        operator()(Args&&... args)
        {
            MARRAY_ASSERT(sizeof...(Args) == dimension());
            pointer ptr = data();
            get_reference<0>(ptr, std::forward<Args>(args)...);
            return *ptr;
        }

        const_pointer cdata() const
        {
            return alloc_;
        }

        const_pointer data() const
        {
            return alloc_;
        }

        pointer data()
        {
            return alloc_;
        }

        idx_type length(unsigned dim) const
        {
            return len_[dim];
        }

        const std::vector<idx_type>& lengths() const
        {
            return len_;
        }

        stride_type stride(unsigned dim) const
        {
            return stride_[dim];
        }

        const std::vector<stride_type>& strides() const
        {
            return stride_;
        }

        unsigned dimension() const
        {
            return static_cast<unsigned>(len_.size());
        }

        void swap(varray& other)
        {
            using std::swap;
            swap(alloc_.data_, other.alloc_.data_);
            swap(alloc_.size_, other.alloc_.size_);
            swap(len_, other.len_);
            swap(stride_, other.stride_);
            swap(layout_, other.layout_);
        }

        friend void swap(varray& a, varray& b)
        {
            a.swap(b);
        }
};

template <unsigned NDim, typename T, typename Alloc>
marray_view<const T, NDim> fix(const varray<T, Alloc>& other)
{
    MARRAY_ASSERT(NDim == other.dimension());

    std::array<idx_type, NDim> len;
    std::array<stride_type, NDim> stride;
    std::copy_n(other.lengths().begin(), NDim, len.begin());
    std::copy_n(other.strides().begin(), NDim, stride.begin());

    return {len, other.data(), stride};
}

template <unsigned NDim, typename T, typename Alloc>
marray_view<T, NDim> fix(varray<T, Alloc>& other)
{
    MARRAY_ASSERT(NDim == other.dimension());

    std::array<idx_type, NDim> len;
    std::array<stride_type, NDim> stride;
    std::copy_n(other.lengths().begin(), NDim, len.begin());
    std::copy_n(other.strides().begin(), NDim, stride.begin());

    return {len, other.data(), stride};
}

template <typename T, typename U>
detail::enable_if_assignable_t<U&,T>
copy(const varray_view<T>& a, const varray_view<U>& b)
{
    MARRAY_ASSERT(a.lengths() == b.lengths());

    auto it = make_iterator(a.lengths(), a.strides(), b.strides());
    auto a_ = a.data();
    auto b_ = b.data();
    while (it.next(a_, b_)) *b_ = *a_;
}

template <typename T, typename U>
detail::enable_if_assignable_t<U&,T>
copy(const T& a, const varray_view<U>& b)
{
    auto it = make_iterator(b.lengths(), b.strides());
    auto b_ = b.data();
    while (it.next(b_)) *b_ = a;
}

}

#endif
