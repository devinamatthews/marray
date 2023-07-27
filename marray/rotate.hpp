#ifndef MARRAY_ROTATE_HPP
#define MARRAY_ROTATE_HPP

#include "marray.hpp"

namespace MArray
{

template <typename T>
std::enable_if_t<detail::is_marray_like_v<T>>
rotate(T&& array, const array_1d<len_type>& shift)
{
    MARRAY_ASSERT(shift.size() == array.dimension());

    for (auto i : range(array.dimension()))
        rotate(array, i, shift[i]);
}

template <typename T>
std::enable_if_t<detail::is_marray_like_v<T>>
rotate(T&& array, int dim, len_type shift)
{
    MARRAY_ASSERT(dim >= 0 && dim < array.dimension());

    auto n = array.length(dim);
    auto s = array.stride(dim);

    if (n == 0) return;

    shift = shift%n;
    if (shift < 0) shift += n;

    if (shift == 0) return;

    auto len = array.lengths();
    auto& stride = array.strides();
    len[dim] = 1;

    auto p = array.data();
    auto it = make_iterator(len, stride);
    while (it.next(p))
    {
        auto a = p;
        auto b = p+(shift-1)*s;
        while (a < b)
        {
            std::iter_swap(a, b);
            a += s;
            b -= s;
        }

        a = p+shift*s;
        b = p+(n-1)*s;
        while (a < b)
        {
            std::iter_swap(a, b);
            a += s;
            b -= s;
        }

        a = p;
        b = p+(n-1)*s;
        while (a < b)
        {
            std::iter_swap(a, b);
            a += s;
            b -= s;
        }
    }
}

}

#endif //MARRAY_ROTATE_HPP
