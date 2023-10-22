#ifndef MARRAY_ARRAY_1D_HPP
#define MARRAY_ARRAY_1D_HPP

#include <type_traits>
#include <initializer_list>
#include <vector>
#include <array>

#include "detail/utility.hpp"

namespace MArray
{
namespace detail
{

template <typename T>
struct is_row : std::false_type {};

template <typename T, typename D, bool O, int Tags>
struct is_row<marray_base<T, 1, D, O, Tags>> : std::true_type {};

template <typename T, int Tags>
struct is_row<marray_view<T, 1, Tags>> : std::true_type {};

template <typename T, typename A>
struct is_row<marray<T, 1, A>> : std::true_type {};

template <typename T, typename U, typename=void>
struct is_row_of : std::false_type {};

template <typename T, typename U>
struct is_row_of<T, U, std::enable_if_t<is_row<T>::value &&
    std::is_assignable<U&,typename T::value_type>::value>> : std::true_type {};

template <typename T, typename U>
struct is_1d_container_of :
    std::integral_constant<bool, is_row_of<T,U>::value ||
                                 is_container_of<T,U>::value> {};

template <typename T, typename U, typename V=void>
using enable_if_1d_container_of_t = std::enable_if_t<is_1d_container_of<T,U>::value, V>;

template <typename T>
std::enable_if_t<is_container<T>::value,len_type>
length(const T& len)
{
    return len.size();
}

template <typename T>
std::enable_if_t<is_row<T>::value,len_type>
length(const T& len)
{
    return len.length();
}

}

/**
 * Adaptor class which can capture a container, initializer_list, @ref row, or @ref row_view with elements
 * convertible to a specified type.
 *
 * This class is used to capture a variety of acceptable argument types without relying on template parameters.
 * For example:
 *
 * @code{.cxx}
 * void foo(const array_1d<int>&);
 *
 * // An initializer_list<int> (sort of...actually narrowing conversions are allowed):
 * foo({1, 4, 4l, 3u});
 *
 * // A linear container:
 * std::vector<long> v{1, 5, 2};
 * foo(v);
 *
 * // A non-linear container:
 * std::list<int> l{1, 6, 2, 0};
 * foo(l);
 *
 * // A row (1-d tensor):
 * row<int> r{3};
 * r[0] = 0;
 * r[1] = 4;
 * r[2] = 2;
 * foo(r);
 *
 * // A row_view:
 * row_view<int> rv({v.size()}, v.data());
 * foo(rv);
 * @endcode
 *
 * @tparam T  Elements of the supplied container, initializer list, etc. must be convertible to this type.
 *
 * @ingroup util
 */
template <typename T, int NDim=DYNAMIC>
class array_1d : public detail::array_type_t<T, NDim>
{
    protected:
        using base = detail::array_type_t<T, NDim>;

        struct len_wrapper
        {
            T len;

            template <typename U, typename=std::enable_if_t<std::is_convertible_v<U,T>>>
            len_wrapper(const U& len) : len{(T)len} {}

            operator T() const { return len; }
        };

    public:
        array_1d() {}

        template <int N=NDim>
        array_1d(std::initializer_list<len_wrapper> il, std::enable_if_t<N == DYNAMIC>* = nullptr)
        : base(il.begin(), il.end()) {}

        template <int N=NDim>
        array_1d(std::initializer_list<len_wrapper> il, std::enable_if_t<N != DYNAMIC>* = nullptr)
        {
            MARRAY_ASSERT(il.size() == NDim);
            std::copy_n(il.begin(), NDim, this->data());
        }

        template <typename U, int N=NDim, typename=std::enable_if_t<detail::is_1d_container_of<U,T>::value &&
                                                                    !std::is_same_v<U,std::initializer_list<len_wrapper>>>>
        array_1d(const U& data, std::enable_if_t<N == DYNAMIC>* = nullptr)
        : base(data.begin(), data.end()) {}

        template <typename U, int N=NDim, typename=std::enable_if_t<detail::is_1d_container_of<U,T>::value &&
                                                                    !std::is_same_v<U,std::initializer_list<len_wrapper>>>>
        array_1d(const U& data, std::enable_if_t<N != DYNAMIC>* = nullptr)
        {
            MARRAY_ASSERT(data.size() == NDim);
            std::copy_n(data.begin(), NDim, this->data());
        }

        void slurp(base& other)
        {
            other = std::move(*this);
        }

        void slurp(base& other) const
        {
            other = *this;
        }
};

}

#endif //MARRAY_ARRAY_1D_HPP
