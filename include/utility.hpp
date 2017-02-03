#ifndef _MARRAY_UTILITY_HPP_
#define _MARRAY_UTILITY_HPP_

#include <type_traits>
#include <array>
#include <vector>
#include <utility>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <cmath>
#include <cstddef>

#ifndef MARRAY_DEFAULT_LAYOUT
#define MARRAY_DEFAULT_LAYOUT ROW_MAJOR
#endif

#ifdef MARRAY_ENABLE_ASSERTS
#define MARRAY_ASSERT(e) assert(e)
#else
#define MARRAY_ASSERT(e)
#endif

namespace MArray
{

/*
 * The type all_t specifies a range [0,len_i) for an array
 * dimension i of length len_i (i.e. it selects all of the data along
 * that dimension).
 */
struct all_t { constexpr all_t() {} };
struct bcast_t { constexpr bcast_t() {} };
namespace slice
{
    constexpr all_t all;
    constexpr bcast_t bcast;
}

typedef ptrdiff_t idx_type;
typedef size_t size_type;
typedef ptrdiff_t stride_type;

/*
 * The special value uninitialized is used to construct an array which
 * does not default- or value-initialize its elements (useful for avoiding
 * redundant memory operations for scalar types).
 */
struct uninitialized_t { constexpr uninitialized_t() {} };
constexpr uninitialized_t uninitialized;

/*
 * Specifies the layout of the array data.
 */
enum class layout {COLUMN_MAJOR, ROW_MAJOR, DEFAULT=MARRAY_DEFAULT_LAYOUT};

template <typename I> class range_t;

namespace detail
{
    template <typename T>
    using decay_t = typename std::decay<T>::type;

    template <typename T>
    using remove_cv_t = typename std::remove_cv<T>::type;

    template <bool Cond, typename T=void>
    using enable_if_t = typename std::enable_if<Cond,T>::type;

    template <bool Cond, typename T, typename U>
    using conditional_t = typename std::conditional<Cond,T,U>::type;

    template <typename T, typename U=void>
    using enable_if_integral_t = enable_if_t<std::is_integral<T>::value,U>;

    template <typename T, typename U=void>
    using enable_if_not_integral_t = enable_if_t<!std::is_integral<T>::value,U>;

    template <typename T, typename U=void>
    using enable_if_const_t = enable_if_t<std::is_const<T>::value,U>;

    template <typename T, typename U=void>
    using enable_if_not_const_t = enable_if_t<!std::is_const<T>::value,U>;

    template <typename T, typename U, typename V=void>
    using enable_if_assignable_t = enable_if_t<std::is_assignable<T,U>::value,V>;

    template <typename T, typename U, typename V=void>
    using enable_if_convertible_t = enable_if_t<std::is_convertible<T,U>::value,V>;

    template <typename T, typename... Args>
    struct are_convertible;

    template <typename T>
    struct are_convertible<T> : std::true_type {};

    template <typename T, typename Arg, typename... Args>
    struct are_convertible<T, Arg, Args...> :
        conditional_t<std::is_convertible<Arg, T>::value,
                      are_convertible<T, Args...>,
                      std::false_type> {};

    template <typename T, typename=void>
    struct is_index_or_slice_helper : std::false_type {};

    template <typename T>
    struct is_index_or_slice_helper<T, enable_if_convertible_t<T, int>> : std::true_type {};

    template <typename I>
    struct is_index_or_slice_helper<range_t<I>, enable_if_integral_t<I>> : std::true_type {};

    template <>
    struct is_index_or_slice_helper<all_t> : std::true_type {};

    template <typename T>
    struct is_index_or_slice : is_index_or_slice_helper<typename std::decay<T>::type> {};

    template <typename... Args>
    struct are_indices_or_slices;

    template<>
    struct are_indices_or_slices<> : std::true_type {};

    template <typename Arg, typename... Args>
    struct are_indices_or_slices<Arg, Args...> :
        conditional_t<is_index_or_slice<Arg>::value,
                      are_indices_or_slices<Args...>,
                      std::false_type> {};

    template <typename...>
    struct is_container_helper {};

    template <typename T, typename=void>
    struct is_container : std::false_type {};

    template <typename T>
    struct is_container<T,
        typename std::conditional<false,
                                  is_container_helper<typename T::value_type,
                                                      decltype(std::declval<T>().size()),
                                                      decltype(std::declval<T>().begin()),
                                                      decltype(std::declval<T>().end())>,
                                  void>::type>
    : std::true_type {};

    template <typename T, typename C, typename=void>
    struct is_container_of : std::false_type {};

    template <typename T, typename C>
    struct is_container_of<T, C, typename std::enable_if<is_container<C>::value>::type>
    : std::is_convertible<typename C::value_type, T> {};

    template <typename T, typename... Ts>
    struct are_containers_helper;

    template <typename T>
    struct are_containers_helper<T> : is_container<T> {};

    template <typename T, typename... Ts>
    struct are_containers_helper
    : std::conditional<is_container<T>::value,
                       are_containers_helper<Ts...>,
                       std::false_type>::type {};

    template <typename... Ts>
    struct are_containers;

    template <>
    struct are_containers<> : std::true_type {};

    template <typename... Ts>
    struct are_containers : are_containers_helper<Ts...> {};

    template <typename T, typename C, typename... Cs>
    struct are_containers_of_helper;

    template <typename T, typename C>
    struct are_containers_of_helper<T, C> : is_container_of<T, C> {};

    template <typename T, typename C, typename... Cs>
    struct are_containers_of_helper
    : std::conditional<is_container_of<T, C>::value,
                       are_containers_of_helper<T, Cs...>,
                       std::false_type>::type {};

    template <typename T, typename... Cs>
    struct are_containers_of;

    template <typename T>
    struct are_containers_of<T> : std::true_type {};

    template <typename T, typename... Cs>
    struct are_containers_of : are_containers_of_helper<T, Cs...> {};

    template <typename T, size_t N, typename C>
    struct is_array_of : std::false_type {};

    template <typename T, size_t N>
    struct is_array_of<T, N, std::array<T,N>> : std::true_type {};

    template <typename T, size_t N, typename C, typename... Cs>
    struct are_arrays_of_helper;

    template <typename T, size_t N, typename C>
    struct are_arrays_of_helper<T, N, C> : is_array_of<T, N, C> {};

    template <typename T, size_t N, typename C, typename... Cs>
    struct are_arrays_of_helper
    : std::conditional<is_array_of<T, N, C>::value,
                       are_arrays_of_helper<T, N, Cs...>,
                       std::false_type>::type {};

    template <typename T, size_t N, typename... Cs>
    struct are_arrays_of;

    template <typename T, size_t N>
    struct are_arrays_of<T, N> : std::true_type {};

    template <typename T, size_t N, typename... Cs>
    struct are_arrays_of : are_arrays_of_helper<T, N, Cs...> {};

    template <typename T, typename C>
    struct is_vector_of : std::false_type {};

    template <typename T>
    struct is_vector_of<T, std::vector<T>> : std::true_type {};

    template <typename T, typename C, typename... Cs>
    struct are_vectors_of_helper;

    template <typename T, typename C>
    struct are_vectors_of_helper<T, C> : is_vector_of<T, C> {};

    template <typename T, typename C, typename... Cs>
    struct are_vectors_of_helper
    : std::conditional<is_vector_of<T, C>::value,
                       are_vectors_of_helper<T, Cs...>,
                       std::false_type>::type {};

    template <typename T, typename... Cs>
    struct are_vectors_of;

    template <typename T>
    struct are_vectors_of<T> : std::true_type {};

    template <typename T, typename... Cs>
    struct are_vectors_of : are_vectors_of_helper<T, Cs...> {};

    template <size_t NDim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct inc_offsets_helper;

    template <size_t NDim, size_t N, typename Offset>
    struct inc_offsets_helper<NDim, N, N, Offset>
    {
        template <typename stride_type>
        inc_offsets_helper(unsigned i,
                           Offset& off0,
                           const std::array<std::array<stride_type,NDim>,N>& strides)
        {
            off0 += strides[N-1][i];
        }

        template <typename stride_type>
        inc_offsets_helper(unsigned i,
                           Offset& off0,
                           const std::array<std::vector<stride_type>,N>& strides)
        {
            off0 += strides[N-1][i];
        }
    };

    template <size_t NDim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct inc_offsets_helper
    {
        template <typename stride_type>
        inc_offsets_helper(unsigned i,
                           Offset& off0, Offsets&... off,
                           const std::array<std::array<stride_type,NDim>,N>& strides)
        {
            off0 += strides[I-1][i];
            inc_offsets_helper<NDim, N, I+1, Offsets...>(i, off..., strides);
        }

        template <typename stride_type>
        inc_offsets_helper(unsigned i,
                           Offset& off0, Offsets&... off,
                           const std::array<std::vector<stride_type>,N>& strides)
        {
            off0 += strides[I-1][i];
            inc_offsets_helper<0, N, I+1, Offsets...>(i, off..., strides);
        }
    };

    template <typename stride_type, size_t NDim, size_t N, typename... Offsets>
    void inc_offsets(unsigned i,
                     const std::array<std::array<stride_type,NDim>,N>& strides,
                     Offsets&... off)
    {
        inc_offsets_helper<NDim, N, 1, Offsets...>(i, off..., strides);
    }

    template <typename stride_type, size_t N, typename... Offsets>
    void inc_offsets(unsigned i,
                     const std::array<std::vector<stride_type>,N>& strides,
                     Offsets&... off)
    {
        inc_offsets_helper<0, N, 1, Offsets...>(i, off..., strides);
    }

    template <size_t NDim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct dec_offsets_helper;

    template <size_t NDim, size_t N, typename Offset>
    struct dec_offsets_helper<NDim, N, N, Offset>
    {
        template <typename idx_type, typename stride_type>
        dec_offsets_helper(unsigned i,
                           Offset& off0,
                           const std::array<idx_type,NDim>& pos,
                           const std::array<std::array<stride_type,NDim>,N>& strides)
        {
            off0 -= pos[i]*strides[N-1][i];
        }

        template <typename idx_type, typename stride_type>
        dec_offsets_helper(unsigned i,
                           Offset& off0,
                           const std::vector<idx_type>& pos,
                           const std::array<std::vector<stride_type>,N>& strides)
        {
            off0 -= pos[i]*strides[N-1][i];
        }
    };

    template <size_t NDim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct dec_offsets_helper
    {
        template <typename idx_type, typename stride_type>
        dec_offsets_helper(unsigned i,
                           Offset& off0, Offsets&... off,
                           const std::array<idx_type,NDim>& pos,
                           const std::array<std::array<stride_type,NDim>,N>& strides)
        {
            off0 -= pos[i]*strides[I-1][i];
            dec_offsets_helper<NDim, N, I+1, Offsets...>(i, off..., pos, strides);
        }

        template <typename idx_type, typename stride_type>
        dec_offsets_helper(unsigned i,
                           Offset& off0, Offsets&... off,
                           const std::vector<idx_type>& pos,
                           const std::array<std::vector<stride_type>,N>& strides)
        {
            off0 -= pos[i]*strides[I-1][i];
            dec_offsets_helper<0, N, I+1, Offsets...>(i, off..., pos, strides);
        }
    };

    template <typename idx_type, typename stride_type, size_t NDim, size_t N, typename... Offsets>
    void dec_offsets(unsigned i,
                     const std::array<idx_type,NDim>& pos,
                     const std::array<std::array<stride_type,NDim>,N>& strides,
                     Offsets&... off)
    {
        dec_offsets_helper<NDim, N, 1, Offsets...>(i, off..., pos, strides);
    }

    template <typename idx_type, typename stride_type, size_t N, typename... Offsets>
    void dec_offsets(unsigned i,
                     const std::vector<idx_type>& pos,
                     const std::array<std::vector<stride_type>,N>& strides,
                     Offsets&... off)
    {
        dec_offsets_helper<0, N, 1, Offsets...>(i, off..., pos, strides);
    }

    template <size_t NDim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct move_offsets_helper;

    template <size_t NDim, size_t N, typename Offset>
    struct move_offsets_helper<NDim, N, N, Offset>
    {
        template <typename idx_type, typename stride_type>
        move_offsets_helper(Offset& off0,
                           const std::array<idx_type,NDim>& pos,
                           const std::array<std::array<stride_type,NDim>,N>& strides)
        {
            for (unsigned i = 0;i < pos.size();i++) off0 += pos[i]*strides[N-1][i];
        }

        template <typename idx_type, typename stride_type>
        move_offsets_helper(Offset& off0,
                           const std::vector<idx_type>& pos,
                           const std::array<std::vector<stride_type>,N>& strides)
        {
            for (unsigned i = 0;i < pos.size();i++) off0 += pos[i]*strides[N-1][i];
        }
    };

    template <size_t NDim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct move_offsets_helper
    {
        template <typename idx_type, typename stride_type>
        move_offsets_helper(Offset& off0, Offsets&... off,
                           const std::array<idx_type,NDim>& pos,
                           const std::array<std::array<stride_type,NDim>,N>& strides)
        {
            for (unsigned i = 0;i < pos.size();i++) off0 += pos[i]*strides[I-1][i];
            move_offsets_helper<NDim, N, I+1, Offsets...>(off..., pos, strides);
        }

        template <typename idx_type, typename stride_type>
        move_offsets_helper(Offset& off0, Offsets&... off,
                           const std::vector<idx_type>& pos,
                           const std::array<std::vector<stride_type>,N>& strides)
        {
            for (unsigned i = 0;i < pos.size();i++) off0 += pos[i]*strides[I-1][i];
            move_offsets_helper<0, N, I+1, Offsets...>(off..., pos, strides);
        }
    };

    template <typename idx_type, typename stride_type, size_t NDim, size_t N, typename... Offsets>
    void move_offsets(const std::array<idx_type,NDim>& pos,
                     const std::array<std::array<stride_type,NDim>,N>& strides,
                     Offsets&... off)
    {
        move_offsets_helper<NDim, N, 1, Offsets...>(off..., pos, strides);
    }

    template <typename idx_type, typename stride_type, size_t N, typename... Offsets>
    void move_offsets(const std::vector<idx_type>& pos,
                     const std::array<std::vector<stride_type>,N>& strides,
                     Offsets&... off)
    {
        move_offsets_helper<0, N, 1, Offsets...>(off..., pos, strides);
    }

    template <size_t NDim, size_t N, size_t I, typename Stride, typename... Strides>
    struct set_strides_helper;

    template <size_t NDim, size_t N, typename Stride>
    struct set_strides_helper<NDim, N, N, Stride>
    {
        template <typename stride_type>
        set_strides_helper(const Stride& stride0,
                           std::array<std::array<stride_type,NDim>,N>& strides_)
        {
            assert(stride0.size() == NDim);
            std::copy_n(stride0.begin(), NDim, strides_[N-1].begin());
        }

        template <typename stride_type>
        set_strides_helper(const Stride& stride0,
                           std::array<std::vector<stride_type>,N>& strides_)
        {
            strides_[N-1].assign(stride0.begin(), stride0.end());
        }
    };

    template <size_t NDim, size_t N, size_t I, typename Stride, typename... Strides>
    struct set_strides_helper
    {
        template <typename stride_type>
        set_strides_helper(const Stride& stride0, const Strides&... strides,
                           std::array<std::array<stride_type,NDim>,N>& strides_)
        {
            assert(stride0.size() == NDim);
            std::copy_n(stride0.begin(), NDim, strides_[I-1].begin());
            set_strides_helper<NDim, N, I+1, Strides...>(strides..., strides_);
        }

        template <typename stride_type>
        set_strides_helper(const Stride& stride0, const Strides&... strides,
                           std::array<std::vector<stride_type>,N>& strides_)
        {
            strides_[I-1].assign(stride0.begin(), stride0.end());
            set_strides_helper<0, N, I+1, Strides...>(strides..., strides_);
        }
    };

    template <typename stride_type, size_t NDim, size_t N, typename... Strides>
    void set_strides(std::array<std::array<stride_type,NDim>,N>& strides_,
                     const Strides&... strides)
    {
        set_strides_helper<NDim, N, 1, Strides...>(strides..., strides_);
    }

    template <typename stride_type, size_t N, typename... Strides>
    void set_strides(std::array<std::vector<stride_type>,N>& strides_,
                     const Strides&... strides)
    {
        set_strides_helper<0, N, 1, Strides...>(strides..., strides_);
    }

    template <typename T, T... S> struct integer_sequence {};

    template <typename T, typename U, typename V> struct concat_sequences;
    template <typename T, T... S, T... R>
    struct concat_sequences<T, integer_sequence<T, S...>, integer_sequence<T, R...>>
    {
        typedef integer_sequence<T, S..., (R+sizeof...(S))...> type;
    };

    template <size_t N> struct static_range_helper;

    template <> struct static_range_helper<0>
    {
        typedef integer_sequence<size_t> type;
    };

    template <> struct static_range_helper<1>
    {
        typedef integer_sequence<size_t,0> type;
    };

    template <size_t N> struct static_range_helper
    {
        typedef typename concat_sequences<size_t, typename static_range_helper<(N+1)/2>::type,
                                                  typename static_range_helper<N/2>::type>::type type;
    };

    template <size_t N>
    using static_range = typename static_range_helper<N>::type;
}

/*
 * Create a vector from the specified elements, where the type of the vector
 * is taken from the first element.
 */
template <typename T, typename... Args>
std::vector<typename std::decay<T>::type>
make_vector(T&& t, Args&&... args)
{
    return {{std::forward<T>(t), std::forward<Args>(args)...}};
}

/*
 * Create an array from the specified elements, where the type of the array
 * is taken from the first element.
 */
template <typename T, typename... Args>
std::array<typename std::decay<T>::type, sizeof...(Args)+1>
make_array(T&& t, Args&&... args)
{
    return {{std::forward<T>(t), std::forward<Args>(args)...}};
}

template <typename... Old, typename New>
std::tuple<Old..., detail::decay_t<New>>
push_back(const std::tuple<Old...>& x, New&& y)
{
    return std::tuple_cat(x, std::make_tuple(std::forward<New>(y)));
}

template <typename... Old, typename New>
std::tuple<Old..., detail::decay_t<New>>
push_back(std::tuple<Old...>&& x, New&& y)
{
    return std::tuple_cat(std::move(x), std::make_tuple(std::forward<New>(y)));
}

}

#endif
