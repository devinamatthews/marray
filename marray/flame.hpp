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

struct direction
{
    int dir;

    constexpr direction(int dir) : dir(dir) {}

    bool operator==(const direction& other) const { return dir == other.dir; }
    bool operator!=(const direction& other) const { return dir != other.dir; }
};

constexpr inline direction DOWN{0};
constexpr inline direction UP{1};

constexpr inline const direction& RIGHT = DOWN;
constexpr inline const direction& LEFT = UP;

constexpr inline const direction& BOTTOM_RIGHT = DOWN;
constexpr inline const direction& TOP_LEFT = UP;

constexpr inline const direction& FORWARD = DOWN;
constexpr inline const direction& BACKWARD = UP;

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

template <typename I>
range_t<len_type> convert(const range_t<I>& x)
{
    return x;
}

inline len_type convert(len_type i) { return i; }

template <typename I>
I front(const range_t<I>& x) { return x.front(); };
template <typename I>
I back(const range_t<I>& x) { return x.back(); };
template <typename I>
I size(const range_t<I>& x) { return x.size(); };

inline len_type front(len_type x) { return x; };
inline len_type back(len_type x) { return x; };
inline len_type size(len_type) { return 1; };

template <typename T> struct is_range : std::false_type {};
template <typename I> struct is_range<range_t<I>> : std::true_type {};

template <typename T> constexpr auto is_range_v = is_range<T>::value;

template <size_t I, typename... Args> using nth_type = std::tuple_element_t<I, std::tuple<Args...>>;

template <size_t I, typename... Args> decltype(auto) nth_arg(Args&&... args)
{
    return std::get<I>(std::forward_as_tuple(std::forward<Args>(args)...));
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
template <len_type... Sizes, size_t... NewIdx, size_t... OldIdx, size_t... FirstIdx, typename... Args>
auto repartition(len_type bs, direction dir,
                 std::integer_sequence<len_type, Sizes...>,
                 std::index_sequence<NewIdx...>,
                 std::index_sequence<OldIdx...>,
                 std::index_sequence<FirstIdx...>,
                 const Args&... args)
{
    constexpr auto NOld = sizeof...(OldIdx);
    constexpr auto NNew = sizeof...(NewIdx);
    constexpr auto NFixed = ((Sizes == DYNAMIC ? 0 : Sizes) + ...);
    constexpr auto Blocked = ((Sizes == DYNAMIC) + ...) == 1;

    static_assert(NNew >= 1, "At least one exposed range is required");
    static_assert(((Sizes == DYNAMIC) + ...) < 2, "At most one of the exposed ranges may be dynamic");
    static_assert(detail::is_range_v<detail::nth_type<0, Args...>> &&
                  detail::is_range_v<detail::nth_type<NOld-1, Args...>>,
                  "The first and last input ranges must be blocked");

    if (Blocked) MARRAY_ASSERT(bs >= NFixed);
    bs -= NFixed;

    if (dir == FORWARD)
    {
        auto& last = detail::nth_arg<NOld-1>(args...);
        MARRAY_ASSERT(last.size() >= bs + NFixed);

        std::array<len_type,NNew+1> sizes{last.front(), (Sizes == DYNAMIC ? bs : Sizes)...};
        ((sizes[NewIdx+1] += sizes[NewIdx]), ...);

        return std::make_tuple(convert(nth_arg<FirstIdx>(args...))...,
                               make_range<Sizes>(sizes[NewIdx], sizes[NewIdx+1])...,
                               range_t<len_type>{last.front()+bs+NFixed, last.back()+1});
    }
    else
    {
        auto& first = detail::nth_arg<NOld-1>(args...);
        MARRAY_ASSERT(first.size() >= bs + NFixed);

        constexpr std::array rev{Sizes...};

        std::array<len_type,NNew+1> sizes{(rev[NNew-1-NewIdx] == DYNAMIC ? bs : rev[NNew-1-NewIdx])..., first.back()+1};
        (..., (sizes[NewIdx] = sizes[NewIdx+1]-sizes[NewIdx]));

        return std::make_tuple(range_t<len_type>{first.front(), first.back()+1-bs-NFixed},
                               make_range<rev[NNew-1-NewIdx]>(sizes[NewIdx], sizes[NewIdx+1])...,
                               convert(nth_arg<FirstIdx+1>(args...))...);
    }
}

template <size_t N, len_type... Sizes, typename... Args>
auto repartition(len_type bs, direction dir, const Args&... args)
{
    constexpr auto NNew = sizeof...(Sizes);
    static_assert(N >= 2, "At least two input ranges are required");

    return repartition(bs, dir,
                       std::integer_sequence<len_type, Sizes...>{},
                       std::make_index_sequence<NNew>{},
                       std::make_index_sequence<N>{},
                       std::make_index_sequence<N-1>{},
                       args...);
}

template <size_t... Idx, size_t... FirstIdx, typename... Args>
auto continue_with(direction dir,
                   std::index_sequence<Idx...>,
                   std::index_sequence<FirstIdx...>,
                   const Args&... args)
{
    constexpr auto NOld = sizeof...(Idx);
    constexpr auto NNew = sizeof...(FirstIdx)+1;

    static_assert(detail::is_range_v<detail::nth_type<0, Args...>> &&
                  detail::is_range_v<detail::nth_type<NOld-1, Args...>>,
                  "The first and last input ranges must be blocked");

    auto first = nth_arg<0>(args...).front();
    auto last = nth_arg<NOld-1>(args...).back();
    (MARRAY_ASSERT(back(nth_arg<Idx == NOld-1 ? NOld-2 : Idx>(args...))+1 ==
                   front(nth_arg<Idx == NOld-1 ? NOld-1 : Idx+1>(args..., last+1))), ...);

    if (dir == FORWARD)
    {
        return std::make_tuple(range(first, back(nth_arg<NOld-NNew>(args...))+1), nth_arg<NOld-NNew+1+FirstIdx>(args...)...);
    }
    else
    {
        return std::make_tuple(nth_arg<FirstIdx>(args...)..., range(front(nth_arg<NNew-1>(args...)), last+1));
    }
}

} // namespace detail

// Partitioning

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

// Row and column partitioning

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

// Blocked/unblocked repartition

template <len_type... Sizes, typename... Args>
auto repartition(const Args&... args)
{
    constexpr auto NNew = sizeof...(Sizes);
    constexpr auto NArgs = sizeof...(args);
    constexpr auto NFixed = ((Sizes == DYNAMIC ? 0 : Sizes) + ... + 0);

    auto bs = NFixed+1;
    auto dir = FORWARD;

    if constexpr (NArgs > 1 && std::is_same_v<detail::nth_type<NArgs-1, Args...>, direction> &&
                               std::is_convertible_v<detail::nth_type<NArgs-2, Args...>, len_type>)
    {
        dir = detail::nth_arg<NArgs-1>(args...);
        bs = detail::nth_arg<NArgs-2>(args...);

        if constexpr (NNew == 0)
            return detail::repartition<NArgs-2, DYNAMIC>(bs, dir, args...);
        else
            return detail::repartition<NArgs-2, Sizes...>(bs, dir, args...);
    }
    else if constexpr (NArgs > 1 && std::is_same_v<detail::nth_type<NArgs-2, Args...>, direction> &&
                                    std::is_convertible_v<detail::nth_type<NArgs-1, Args...>, len_type>)
    {
        dir = detail::nth_arg<NArgs-2>(args...);
        bs = detail::nth_arg<NArgs-1>(args...);

        if constexpr (NNew == 0)
            return detail::repartition<NArgs-2, DYNAMIC>(bs, dir, args...);
        else
            return detail::repartition<NArgs-2, Sizes...>(bs, dir, args...);
    }
    else if constexpr (NArgs > 0 && std::is_same_v<detail::nth_type<NArgs-1, Args...>, direction>)
    {
        dir = detail::nth_arg<NArgs-1>(args...);

        if constexpr (NNew == 0)
            return detail::repartition<NArgs-1, 1>(bs, dir, args...);
        else
            return detail::repartition<NArgs-1, Sizes...>(bs, dir, args...);
    }
    else if constexpr (NArgs > 0 && std::is_convertible_v<detail::nth_type<NArgs-1, Args...>, len_type>)
    {
        bs = detail::nth_arg<NArgs-1>(args...);

        if constexpr (NNew == 0)
            return detail::repartition<NArgs-1, DYNAMIC>(bs, dir, args...);
        else
            return detail::repartition<NArgs-1, Sizes...>(bs, dir, args...);
    }
    else
    {
        if constexpr (NNew == 0)
            return detail::repartition<NArgs, 1>(bs, dir, args...);
        else
            return detail::repartition<NArgs, Sizes...>(bs, dir, args...);
    }
}

// Blocked/unblocked continue with

template <size_t N, typename... Args>
auto continue_with(const Args&... args)
{
    constexpr auto NArgs = sizeof...(args);

    static_assert(N >= 1, "At least one range must be merged");

    if constexpr (NArgs > 0 && std::is_same_v<detail::nth_type<NArgs-1, Args...>, direction>)
    {
        static_assert(NArgs >= N+3, "At least two ranges must left after merging");
        return detail::continue_with(detail::nth_arg<NArgs-1>(args...),
                                     std::make_index_sequence<NArgs-1>{},
                                     std::make_index_sequence<NArgs-N-2>{},
                                     args...);
    }
    else
    {
        static_assert(NArgs >= N+2, "At least two ranges must left after merging");
        return detail::continue_with(FORWARD,
                                     std::make_index_sequence<NArgs>{},
                                     std::make_index_sequence<NArgs-N-1>{},
                                     args...);
    }
}

template <typename... Args>
auto continue_with(const Args&... args)
{
    return continue_with<1>(args...);
}

// Diagonal extraction

template <typename MArray>
auto diag(const MArray& A, len_type off=0)
{
    MARRAY_ASSERT(A.dimension() == 2);

    using T = typename MArray::value_type;

    auto m = A.length(0);
    auto n = A.length(1);

    auto m_begin = std::max(len_type{}, off);
    auto n_begin = std::max(len_type{}, -off);

    MARRAY_ASSERT(m_begin <= m);
    MARRAY_ASSERT(n_begin <= n);

    return marray_view<T>{{std::min(m - m_begin, n - n_begin)},
                          A.data() + m_begin*A.stride(0) + n_begin*A.stride(1),
                          {A.stride(0) + A.stride(1)}};
}

template <typename MArray>
auto subdiag(const MArray& A)
{
    return diag(A, 1);
}

template <typename MArray>
void pivot_rows(const MArray& A, len_type pi)
{
    MARRAY_ASSERT(A.dimension() == 2);
    MARRAY_ASSERT(pi >= 0 && pi < A.length(0));

    if (pi == 0) return;

    blas::swap(A.length(1), A.data(), A.stride(1),
               A.data() + pi*A.stride(0), A.stride(1));
}

template <typename MArray>
void pivot_columns(const MArray& A, len_type pi)
{
    MARRAY_ASSERT(A.dimension() == 2);
    MARRAY_ASSERT(pi >= 0 && pi < A.length(1));

    if (pi == 0) return;

    blas::swap(A.length(0), A.data(), A.stride(0),
               A.data() + pi*A.stride(1), A.stride(0));
}

template <typename MArray>
void pivot_both(const MArray& A, len_type pi)
{
    pivot_rows(A, pi);
    pivot_colums(A, pi);
}

} //namespace MArray

#endif //MARRAY_FLAME_HPP
