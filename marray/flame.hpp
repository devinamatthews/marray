#ifndef MARRAY_FLAME_HPP
#define MARRAY_FLAME_HPP

#include <type_traits>
#include <utility>
#include <tuple>
#include <array>

#include "marray_view.hpp"
#include "expression.hpp"
#include "blas.h"

namespace MArray
{

enum direction
{
    DOWN,
    UP,
    RIGHT = DOWN,
    LEFT = UP,
    BOTTOM_RIGHT = DOWN,
    UPPER_LEFT = UP,
    FORWARD = DOWN,
    BACKWARD = UP,
};

namespace detail
{

template <typename I, I N>
std::enable_if_t<N == I{1},I> make_range(I from, I)
{
    return from;
}

template <typename I, I N>
std::enable_if_t<N != I{1},range_t<I>> make_range(I from, I to)
{
    return range(from, to);
}

template <typename I, I... Sizes, size_t... Idx>
auto partition(const range_t<I>& x, direction dir, std::integer_sequence<I, Sizes...>, std::index_sequence<Idx...>)
{
    constexpr auto npart = sizeof...(Sizes);
    constexpr auto nfixed = ((Sizes > 0) + ...);
    static_assert(npart >= 2, "At least two partitions are required.");
    static_assert(npart - nfixed == 2, "Exactly two dynamic partitions are required.");

    constexpr std::array sizes0{Sizes...};

    constexpr auto front = [&]
    {
        auto i = 0;
        while (sizes0[i] > 0) i++;
        return i;
    }();

    constexpr auto back = [&]
    {
        auto i = npart-1;
        while (sizes0[i] > 0) i--;
        return i;
    }();

    const auto nleft = x.size() - nfixed;
    std::array sizes{x.front(), Sizes...};
    sizes[front+1] = dir == FORWARD ? 0 : nleft;
    sizes[back+1] = dir == FORWARD ? nleft : 0;

    ((sizes[Idx+1] += sizes[Idx]), ...);

    return std::make_tuple(make_range<I,Sizes>(sizes[Idx], sizes[Idx+1])...);
}

} //namespace detail

//partitioning

template <len_type Size0, len_type... Sizes, typename I>
auto partition(const range_t<I>& x, direction dir = FORWARD)
{
    return detail::partition(x, dir,
                             std::integer_sequence<len_type, Size0, Sizes...>{},
                             std::make_index_sequence<sizeof...(Sizes)+1>{});
}

template <typename I>
auto partition(const range_t<I>& x, direction dir = FORWARD)
{
    return partition<-1, -1>(x, dir);
}

template <len_type... Sizes, typename MArray>
auto partition(const MArray& A, int dim, direction dir = FORWARD)
{
    return partition<Sizes...>(range(A.length(dim)), dir);
}

//row and column partitioning

template <typename MArray>
auto rows(const MArray& A)
{
    return range(A.length(0));
}

template <typename MArray>
auto columns(const MArray& A)
{
    return range(A.length(1));
}

template <len_type... Sizes, typename MArray>
auto partition_rows(const MArray& A, direction dir = FORWARD)
{
    return partition<Sizes...>(A, 0, dir);
}

template <len_type... Sizes, typename MArray>
auto partition_columns(const MArray& A, direction dir = FORWARD)
{
    return partition<Sizes...>(A, 1, dir);
}

// Unblocked repartition

template <typename I>
auto repartition(const range_t<I>& A, const range_t<I>& B, direction dir = FORWARD)
{
    MARRAY_ASSERT(A.back()+1 == B.front());

    if (dir == FORWARD)
    {
        MARRAY_ASSERT(B.size() > 0);
        return std::make_tuple(A, B.from(), range(B.from()+1, B.to()));
    }
    else
    {
        MARRAY_ASSERT(A.size() > 0);
        return std::make_tuple(range(A.from(), A.to()-1), A.to()-1, B);
    }
}

template <typename I>
auto repartition(const range_t<I>& A, I b, const range_t<I>& C, direction dir = FORWARD)
{
    MARRAY_ASSERT(A.back()+1 == b);
    MARRAY_ASSERT(b+1 == C.front());

    if (dir == FORWARD)
    {
        MARRAY_ASSERT(C.size() > 0);
        return std::make_tuple(A, b, b+1, range(C.from()+1, C.to()));
    }
    else
    {
        MARRAY_ASSERT(A.size() > 0);
        return std::make_tuple(range(A.from(), A.to()-1), b-1, b, C);
    }
}

// Blocked repartition

template <typename I>
auto repartition(const range_t<I>& A, const range_t<I>& B, I n, direction dir = FORWARD)
{
    MARRAY_ASSERT(A.back()+1 == B.front());

    if (dir == FORWARD)
    {
        if (n < B.size()) n = B.size();
        return std::make_tuple(A, range(B.front(), B.front()+n), range(B.front()+n, B.back()+1));
    }
    else
    {
        if (n < A.size()) n = A.size();
        return std::make_tuple(range(A.front(), A.back()+1-n), range(A.back()+1-n, A.back()+1), B);
    }
}

// Unblocked continue with

template <typename I>
auto continue_with(const range_t<I>& R0, I r1, const range_t<I>& R2, direction dir = FORWARD)
{
    MARRAY_ASSERT(R0.back()+1 == r1);
    MARRAY_ASSERT(r1+1 == R2.front());

    if (dir == FORWARD)
    {
        return std::make_tuple(range(R0.from(), R0.to()+1), R2);
    }
    else
    {
        return std::make_tuple(R0, range(R2.from()-1, R2.to()));
    }
}

template <typename I>
auto continue_with(const range_t<I>& R0, I r1, I r2, const range_t<I>& R3, direction dir = FORWARD)
{
    MARRAY_ASSERT(R0.back()+1 == r1);
    MARRAY_ASSERT(r1+1 == r2);
    MARRAY_ASSERT(r2+1 == R3.front());

    if (dir == FORWARD)
    {
        return std::make_tuple(range(R0.from(), R0.to()+1), r2, R3);
    }
    else
    {
        return std::make_tuple(R0, r1, range(R3.from()-1, R3.to()));
    }
}

// Blocked continue with

template <typename I>
auto continue_with(const range_t<I>& R0, const range_t<I>& R1, const range_t<I>& R2, direction dir = FORWARD)
{
    MARRAY_ASSERT(R0.back()+1 == R1.front());
    MARRAY_ASSERT(R1.back()+1 == R2.front());

    if (dir == FORWARD)
    {
        return std::make_tuple(range(R0.from(), R1.to()), R2);
    }
    else
    {
        return std::make_tuple(R0, range(R1.from(), R2.to()));
    }
}

} //namespace MArray

#endif //MARRAY_FLAME_HPP
