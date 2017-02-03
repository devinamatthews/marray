#ifndef _MARRAY_ROTATE_HPP_
#define _MARRAY_ROTATE_HPP_

#include "utility.hpp"
#include "miterator.hpp"
#include "viterator.hpp"

namespace MArray
{

template <typename Array>
void rotate_dim(Array& array, unsigned dim, idx_type shift)
{
    MARRAY_ASSERT(dim < array.dimension());

    idx_type n = array.length(dim);
    stride_type s = array.stride(dim);

    if (n == 0) return;

    shift = shift%n;
    if (shift < 0) shift += n;

    if (shift == 0) return;

    auto len = array.lengths();
    auto stride = array.strides();
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

#endif
