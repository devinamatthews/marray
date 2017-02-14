#ifndef _MARRAY_EXPRESSION_HPP_
#define _MARRAY_EXPRESSION_HPP_

#include <x86intrin.h>

#include "utility.hpp"
#include "vector.hpp"

namespace MArray
{

struct bcast_dim {};

struct slice_dim
{
    idx_type len;
    stride_type stride;

    slice_dim(idx_type len, stride_type stride)
    : len(len), stride(stride) {}
};

}

#include "marray.hpp"

namespace MArray
{

template <typename T, typename... Dims>
struct array_expr
{
    typedef T& result_type;

    T* data;
    std::tuple<Dims...> dims;

    template <size_t I>
    using dim_type = typename std::tuple_element<I, std::tuple<Dims...>>::type;

    array_expr(T* data, const Dims&... dims)
    : data(data), dims{dims...} {}

    result_type eval() const
    {
        return *data;
    }

    template <unsigned NDim, unsigned Dim>
    detail::enable_if_t<(Dim < NDim-sizeof...(Dims)),result_type>
    eval_at(idx_type) const
    {
        return *data;
    }

    template <unsigned NDim, unsigned Dim>
    detail::enable_if_t<(Dim >= NDim-sizeof...(Dims)),result_type>
    eval_at(idx_type i) const
    {
        return eval_at(i, std::get<Dim-(NDim-sizeof...(Dims))>(dims));
    }

    result_type eval_at(idx_type i, const slice_dim&) const
    {
        return data[i];
    }

    result_type eval_at(idx_type, const bcast_dim&) const
    {
        return *data;
    }
};

template <typename T>
struct constant_expr
{
    typedef T result_type;

    T data;

    constant_expr(T data) : data(data) {}

    result_type eval() const
    {
        return data;
    }

    template <unsigned NDim, unsigned Dim>
    result_type eval_at(idx_type) const
    {
        return data;
    }
};

template <typename Expr, typename=void> struct is_array_expression;

template <typename Expr, typename=void> struct is_expression;

template <typename Expr, typename=void> struct is_unary_expression;

template <typename Expr, typename=void> struct is_binary_expression;

template <typename Expr, typename=void> struct expr_result_type;

template <typename Expr>
struct expr_result_type<Expr, detail::enable_if_t<std::is_arithmetic<Expr>::value>>
{
    typedef detail::decay_t<Expr> type;
};

template <typename Expr>
struct expr_result_type<Expr, detail::enable_if_t<is_expression<Expr>::value>>
{
    typedef typename detail::decay_t<Expr>::result_type type;
};

template <typename Expr>
detail::enable_if_t<is_expression<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval(Expr&& expr)
{
    return expr.eval();
}

template <typename Expr>
detail::enable_if_t<std::is_arithmetic<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval(Expr&& expr)
{
    return expr;
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<is_expression<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval_at(Expr&& expr, idx_type i)
{
    return expr.template eval_at<NDim, Dim>(i);
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<std::is_arithmetic<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval_at(Expr&& expr, idx_type i)
{
    return expr;
}

namespace operators
{

struct plus
{
    template <typename T, typename U>
    auto operator()(const T& a, const U& b) const -> decltype(a+b)
    {
        return a+b;
    }
};

struct minus
{
    template <typename T, typename U>
    auto operator()(const T& a, const U& b) const -> decltype(a-b)
    {
        return a-b;
    }
};

struct multiplies
{
    template <typename T, typename U>
    auto operator()(const T& a, const U& b) const -> decltype(a*b)
    {
        return a*b;
    }
};

struct divides
{
    template <typename T, typename U>
    auto operator()(const T& a, const U& b) const -> decltype(a/b)
    {
        return a/b;
    }
};

struct pow
{
    template <typename T, typename U>
    auto operator()(const T& a, const U& b) const -> decltype(std::pow(a,b))
    {
        return std::pow(a,b);
    }
};

struct negate
{
    template <typename T>
    auto operator()(const T& a) const -> decltype(-a)
    {
        return -a;
    }
};

struct exp
{
    template <typename T>
    auto operator()(const T& a) const -> decltype(std::exp(a))
    {
        return std::exp(a);
    }
};

struct sqrt
{
    template <typename T>
    auto operator()(const T& a) const -> decltype(std::sqrt(a))
    {
        return std::sqrt(a);
    }
};

}

template <typename LHS, typename RHS, typename Op>
struct binary_expr
{
    typedef LHS first_type;
    typedef RHS second_type;
    typedef decltype(std::declval<Op>()(
        std::declval<typename expr_result_type<LHS>::type>(),
        std::declval<typename expr_result_type<RHS>::type>())) result_type;

    LHS first;
    RHS second;
    Op op;

    binary_expr(const LHS& first, const RHS& second, const Op& op = Op())
    : first(first), second(second), op(op) {}

    result_type eval() const
    {
        return op(MArray::eval(first), MArray::eval(second));
    }

    template <unsigned NDim, unsigned Dim>
    result_type eval_at(idx_type i) const
    {
        return op(MArray::eval_at<NDim, Dim>(first, i),
                  MArray::eval_at<NDim, Dim>(second, i));
    }
};

template <typename LHS, typename RHS>
using add_expr = binary_expr<LHS, RHS, operators::plus>;

template <typename LHS, typename RHS>
using sub_expr = binary_expr<LHS, RHS, operators::minus>;

template <typename LHS, typename RHS>
using mul_expr = binary_expr<LHS, RHS, operators::multiplies>;

template <typename LHS, typename RHS>
using div_expr = binary_expr<LHS, RHS, operators::divides>;

template <typename Base, typename Exponent>
using pow_expr = binary_expr<Base, Exponent, operators::pow>;

template <typename Expr, typename Op>
struct unary_expr
{
    typedef Expr expr_type;
    typedef decltype(std::declval<Op>()(
        std::declval<typename expr_result_type<Expr>::type>())) result_type;

    Expr expr;
    Op op;

    unary_expr(const Expr& expr, const Op& op = Op())
    : expr(expr), op(op) {}

    result_type eval() const
    {
        return op(MArray::eval(expr));
    }

    template <unsigned NDim, unsigned Dim>
    result_type eval_at(idx_type i) const
    {
        return op(MArray::eval_at<NDim, Dim>(expr, i));
    }
};

template <typename Expr>
using negate_expr = unary_expr<Expr, operators::negate>;

template <typename Expr>
using exp_expr = unary_expr<Expr, operators::exp>;

template <typename Expr>
using sqrt_expr = unary_expr<Expr, operators::sqrt>;

template <typename T, typename>
struct is_array_expression : std::false_type {};

template <typename T, typename... Dims>
struct is_array_expression<array_expr<T, Dims...>> : std::true_type {};

template <typename T, typename>
struct is_expression : std::false_type {};

template <typename T, typename... Dims>
struct is_expression<array_expr<T, Dims...>> : std::true_type {};

template <typename T>
struct is_expression<constant_expr<T>> : std::true_type {};

template <typename LHS, typename RHS, typename Op>
struct is_expression<binary_expr<LHS, RHS, Op>> : std::true_type {};

template <typename Expr, typename Op>
struct is_expression<unary_expr<Expr, Op>> : std::true_type {};

template <typename Expr, typename>
struct is_unary_expression : std::false_type {};

template <typename Expr>
struct is_unary_expression<Expr,
    detail::enable_if_t<is_expression<typename Expr::expr_type>::value>>
    : std::true_type {};

template <typename Expr, typename>
struct is_binary_expression : std::false_type {};

template <typename Expr>
struct is_binary_expression<Expr,
    detail::enable_if_t<is_expression<typename Expr::first_type>::value ||
                        is_expression<typename Expr::second_type>::value>>
    : std::true_type {};

template <typename Array, typename Dim>
struct array_expr_type_helper2;

template <typename T, typename... Dims, typename Dim>
struct array_expr_type_helper2<array_expr<T, Dims...>, Dim>
{
    typedef array_expr<T, Dims..., Dim> type;
};

template <typename T, unsigned NDim>
struct array_expr_helper;

template <typename T>
struct array_expr_helper<T, 0>
{
    typedef array_expr<T> type;
};

template <typename T, unsigned NDim>
struct array_expr_helper
{
    typedef typename array_expr_type_helper2<
        typename array_expr_helper<T, NDim-1>::type, slice_dim>::type type;
};

template <typename Array, typename... Dims1>
struct slice_array_expr_helper;

template <typename T, typename... Dims2, typename... Dims1>
struct slice_array_expr_helper<array_expr<T, Dims2...>, Dims1...>
{
    typedef array_expr<T, Dims1..., Dims2...> type;
};

template <typename Expr, typename=void>
struct expression_type;

template <typename T, unsigned NDim, unsigned NIndexed, typename... Dims>
struct expression_type<const marray_slice<T, NDim, NIndexed, Dims...>>
{
    typedef typename slice_array_expr_helper<
        typename array_expr_helper<T, NDim-NIndexed>::type, Dims...>::type type;
};

template <typename T, unsigned NDim, unsigned NIndexed, typename... Dims>
struct expression_type<marray_slice<T, NDim, NIndexed, Dims...>>
{
    typedef typename slice_array_expr_helper<
        typename array_expr_helper<T, NDim-NIndexed>::type, Dims...>::type type;
};

template <typename T, unsigned NDim>
struct expression_type<const marray_view<T, NDim>>
{
    typedef typename array_expr_helper<T, NDim>::type type;
};

template <typename T, unsigned NDim>
struct expression_type<marray_view<T, NDim>>
{
    typedef typename array_expr_helper<T, NDim>::type type;
};

template <typename T, unsigned NDim, typename Alloc>
struct expression_type<const marray<T, NDim, Alloc>>
{
    typedef typename array_expr_helper<const T, NDim>::type type;
};

template <typename T, unsigned NDim, typename Alloc>
struct expression_type<marray<T, NDim, Alloc>>
{
    typedef typename array_expr_helper<T, NDim>::type type;
};

template <typename Expr>
struct expression_type<Expr, detail::enable_if_t<is_expression<typename std::remove_cv<Expr>::type>::value ||
                                                 std::is_arithmetic<Expr>::value>>
{
    typedef typename std::remove_cv<Expr>::type type;
};

template <typename T, unsigned NDim, unsigned NIndexed, typename... Dims,
          size_t... I, size_t... J>
typename expression_type<marray_slice<T, NDim, NIndexed, Dims...>>::type
make_expression_helper(const marray_slice<T, NDim, NIndexed, Dims...>& x,
                       detail::integer_sequence<size_t, I...>,
                       detail::integer_sequence<size_t, J...>)
{
    return {x.data(), x.template dim<I>()...,
            slice_dim(x.template base_length<NIndexed+J>(),
                      x.template base_stride<NIndexed+J>())...};
}

template <typename T, unsigned NDim, size_t... I>
typename expression_type<marray_view<T, NDim>>::type
make_expression_helper(const marray_view<T, NDim>& x,
                       detail::integer_sequence<size_t, I...>)
{
    return {x.data(), slice_dim(x.template length<I>(), x.template stride<I>())...};
}

template <typename T, unsigned NDim, typename Alloc, size_t... I>
typename expression_type<marray<const T, NDim, Alloc>>::type
make_expression_helper(const marray<T, NDim, Alloc>& x,
                       detail::integer_sequence<size_t, I...>)
{
    return {x.data(), slice_dim(x.template length<I>(), x.template stride<I>())...};
}

template <typename T, unsigned NDim, typename Alloc, size_t... I>
typename expression_type<marray<T, NDim, Alloc>>::type
make_expression_helper(marray<T, NDim, Alloc>& x,
                       detail::integer_sequence<size_t, I...>)
{
    return {x.data(), slice_dim(x.template length<I>(), x.template stride<I>())...};
}

template <typename T, unsigned NDim, unsigned NIndexed, typename... Dims>
typename expression_type<marray_slice<T, NDim, NIndexed, Dims...>>::type
make_expression(const marray_slice<T, NDim, NIndexed, Dims...>& x)
{
    return make_expression_helper(x, detail::static_range<sizeof...(Dims)>(),
                                  detail::static_range<NDim-NIndexed>());
}

template <typename T, unsigned NDim>
typename expression_type<marray_view<T, NDim>>::type
make_expression(const marray_view<T, NDim>& x)
{
    return make_expression_helper(x, detail::static_range<NDim>());
}

template <typename T, unsigned NDim, typename Alloc>
typename expression_type<marray<const T, NDim, Alloc>>::type
make_expression(const marray<T, NDim, Alloc>& x)
{
    return make_expression_helper(x, detail::static_range<NDim>());
}

template <typename T, unsigned NDim, typename Alloc>
typename expression_type<marray<T, NDim, Alloc>>::type
make_expression(marray<T, NDim, Alloc>& x)
{
    return make_expression_helper(x, detail::static_range<NDim>());
}

template <typename Expr>
detail::enable_if_t<is_expression<detail::decay_t<Expr>>::value ||
                    std::is_arithmetic<detail::decay_t<Expr>>::value, Expr&&>
make_expression(Expr&& x)
{
    return std::forward<Expr>(x);
}

template <typename Array>
struct is_marray : std::false_type {};

template <typename T, unsigned NDim, unsigned NIndexed, typename... Dims>
struct is_marray<marray_slice<T, NDim, NIndexed, Dims...>> : std::true_type {};

template <typename T, unsigned NDim>
struct is_marray<marray_view<T, NDim>> : std::true_type {};

template <typename T, unsigned NDim, typename Alloc>
struct is_marray<marray<T, NDim, Alloc>> : std::true_type {};

template <typename Expr>
struct is_expression_arg :
    std::integral_constant<bool, is_expression<Expr>::value ||
                                 is_marray<Expr>::value> {};

template <typename Expr>
struct is_expression_arg_or_scalar :
    std::integral_constant<bool, is_expression_arg<Expr>::value ||
                                 std::is_arithmetic<Expr>::value> {};

template <typename LHS, typename RHS>
struct are_expression_args :
    std::integral_constant<bool, is_expression_arg_or_scalar<LHS>::value  &&
                                 is_expression_arg_or_scalar<RHS>::value &&
                                 (!std::is_arithmetic<LHS>::value ||
                                  !std::is_arithmetic<RHS>::value)> {};

template <typename LHS, typename RHS>
detail::enable_if_t<are_expression_args<LHS,RHS>::value,
                    add_expr<typename expression_type<const LHS>::type,
                             typename expression_type<const RHS>::type>>
operator+(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename LHS, typename RHS>
detail::enable_if_t<are_expression_args<LHS,RHS>::value,
                    sub_expr<typename expression_type<const LHS>::type,
                             typename expression_type<const RHS>::type>>
operator-(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename LHS, typename RHS>
detail::enable_if_t<are_expression_args<LHS,RHS>::value,
                    mul_expr<typename expression_type<const LHS>::type,
                             typename expression_type<const RHS>::type>>
operator*(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename LHS, typename RHS>
detail::enable_if_t<are_expression_args<LHS,RHS>::value,
                    div_expr<typename expression_type<const LHS>::type,
                             typename expression_type<const RHS>::type>>
operator/(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename Base, typename Exponent>
detail::enable_if_t<are_expression_args<Base,Exponent>::value,
                    pow_expr<typename expression_type<const Base>::type,
                             typename expression_type<const Exponent>::type>>
pow(const Base& base, const Exponent& exponent)
{
    return {make_expression(base), make_expression(exponent)};
}

template <typename Expr>
detail::enable_if_t<is_expression_arg<Expr>::value,
                    negate_expr<typename expression_type<const Expr>::type>>
operator-(const Expr& expr)
{
    return {make_expression(expr)};
}

template <typename Expr>
detail::enable_if_t<is_expression_arg<Expr>::value,
                    exp_expr<typename expression_type<const Expr>::type>>
exp(const Expr& expr)
{
    return {make_expression(expr)};
}

template <typename Expr>
detail::enable_if_t<is_expression_arg<Expr>::value,
                    sqrt_expr<typename expression_type<const Expr>::type>>
sqrt(const Expr& expr)
{
    return {make_expression(expr)};
}

template <typename Expr, typename=void> struct expr_dimension;

template <typename T, typename... Dims>
struct expr_dimension<array_expr<T, Dims...>>
    : std::integral_constant<unsigned, sizeof...(Dims)> {};

template <typename Expr>
struct expr_dimension<Expr, detail::enable_if_t<std::is_arithmetic<Expr>::value>>
    : std::integral_constant<unsigned, 0> {};

template <typename Expr>
struct expr_dimension<Expr, detail::enable_if_t<is_binary_expression<Expr>::value>>
    : detail::conditional_t<(expr_dimension<typename Expr::first_type>::value <
                             expr_dimension<typename Expr::second_type>::value),
                             expr_dimension<typename Expr::second_type>,
                             expr_dimension<typename Expr::first_type>> {};

template <typename Expr>
struct expr_dimension<Expr, detail::enable_if_t<is_unary_expression<Expr>::value>>
    : expr_dimension<typename Expr::expr_type> {};

template <typename Dim>
detail::enable_if_t<std::is_same<Dim,bcast_dim>::value,idx_type>
get_array_length(const Dim&)
{
    static_assert(!std::is_same<Dim,bcast_dim>::value,
                  "Broadcast dimensions cannot be written to");
    return 0;
}

template <typename Dim>
detail::enable_if_t<std::is_same<Dim,slice_dim>::value,idx_type>
get_array_length(const Dim& dim)
{
    return dim.len;
}

template <typename T, typename... Dims, size_t... I>
std::array<idx_type, sizeof...(I)>
get_array_lengths_helper(const array_expr<T, Dims...>& array,
                         detail::integer_sequence<size_t, I...>)
{
    return {get_array_length(std::get<I>(array.dims))...};
}

template <typename T, typename... Dims>
std::array<idx_type, sizeof...(Dims)>
get_array_lengths(const array_expr<T, Dims...>& array)
{
    return get_array_lengths_helper(array, detail::static_range<sizeof...(Dims)>());
}

inline bool check_expr_length(const bcast_dim&, idx_type)
{
    return true;
}

inline bool check_expr_length(const slice_dim& dim, idx_type len)
{
    return len == dim.len;
}

template <typename T, typename... Dims, size_t NDim, size_t... I>
bool check_expr_lengths_helper(const array_expr<T, Dims...>& array,
                               const std::array<idx_type, NDim>& len,
                               detail::integer_sequence<size_t, I...>)
{
    std::array<bool, sizeof...(Dims)> values =
        {check_expr_length(std::get<I>(array.dims),
                           len[NDim-sizeof...(Dims)+I])...};

    bool ret = true;
    for (bool value : values) ret = ret && value;
    return ret;
}

template <typename T, typename... Dims, size_t NDim>
bool check_expr_lengths(const array_expr<T, Dims...>& array,
                        const std::array<idx_type, NDim>& len)
{
    return check_expr_lengths_helper(array, len, detail::static_range<sizeof...(Dims)>());
}

template <typename Expr, size_t NDim>
detail::enable_if_t<std::is_arithmetic<Expr>::value,bool>
check_expr_lengths(const Expr& expr, const std::array<idx_type, NDim>& len)
{
    return true;
}

template <typename Expr, size_t NDim>
detail::enable_if_t<is_binary_expression<Expr>::value,bool>
check_expr_lengths(const Expr& expr, const std::array<idx_type, NDim>& len)
{
    return check_expr_lengths(expr.first, len) &&
           check_expr_lengths(expr.second, len);
}

template <typename Expr, size_t NDim>
detail::enable_if_t<is_unary_expression<Expr>::value,bool>
check_expr_lengths(const Expr& expr, const std::array<idx_type, NDim>& len)
{
    return check_expr_lengths(expr.expr, len);
}

/*
 * Return true if the dimension is vectorizable (stride-1).
 * Broadcast dimensions are trivially vectorizable.
 */
inline bool is_vectorizable(const bcast_dim& dim)
{
    return true;
}

inline bool is_vectorizable(const slice_dim& dim)
{
    return dim.stride == 1;
}

/*
 * Dim is one of the NDim dimensions of the array being assigned to. Since
 * the number of dimension in the subexpression (sizeof...(Dims)) can be
 * smaller, the initial implicit broadcast dimensions are trivially
 * vectorizable.
 */
template <unsigned NDim, unsigned Dim, typename T, typename... Dims>
detail::enable_if_t<(Dim < NDim-sizeof...(Dims)), bool>
is_vectorizable(array_expr<T, Dims...>& expr)
{
    return true;
}

/*
 * For the remaining sizeof...(Dims) dimensions, subtract NDim-sizeof...(Dims)
 * to get the proper dimension number in the subexpression.
 */
template <unsigned NDim, unsigned Dim, typename T, typename... Dims>
detail::enable_if_t<(Dim >= NDim-sizeof...(Dims)), bool>
is_vectorizable(array_expr<T, Dims...>& expr)
{
    return is_vectorizable(std::get<Dim-(NDim-sizeof...(Dims))>(expr.dims));
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<std::is_arithmetic<Expr>::value, bool>
is_vectorizable(const Expr& expr)
{
    return true;
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<is_binary_expression<Expr>::value, bool>
is_vectorizable(Expr& expr)
{
    return is_vectorizable<NDim, Dim>(expr.first) &&
           is_vectorizable<NDim, Dim>(expr.second);
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<is_unary_expression<Expr>::value, bool>
is_vectorizable(Expr& expr)
{
    return is_vectorizable<NDim, Dim>(expr.expr);
}

/*
 * Increment (go to the next element) or decrement (return to the first element)
 * in a given dimension of the array subexpression.
 *
 * For broadcast dimensions this is a no-op.
 */
template <typename T, typename... Dims>
void increment(array_expr<T, Dims...>& expr, const bcast_dim& dim) {}

template <typename T, typename... Dims>
void increment(array_expr<T, Dims...>& expr, const slice_dim& dim)
{
    expr.data += dim.stride;
}

template <typename T, typename... Dims>
void decrement(array_expr<T, Dims...>& expr, const bcast_dim& dim) {}

template <typename T, typename... Dims>
void decrement(array_expr<T, Dims...>& expr, const slice_dim& dim)
{
    expr.data -= dim.len*dim.stride;
}

/*
 * Dim is one of the NDim dimensions of the array being assigned to. Since
 * the number of dimension in the subexpression (sizeof...(Dims)) can be
 * smaller, ignore increments and decrements for the first NDim-sizeof...(Dims)
 * dimensions.
 */
template <unsigned NDim, unsigned Dim, typename T, typename... Dims>
detail::enable_if_t<(Dim < NDim-sizeof...(Dims))>
increment(array_expr<T, Dims...>& expr) {}

template <unsigned NDim, unsigned Dim, typename T, typename... Dims>
detail::enable_if_t<(Dim < NDim-sizeof...(Dims))>
decrement(array_expr<T, Dims...>& expr) {}

/*
 * For the remaining sizeof...(Dims) dimensions, subtract NDim-sizeof...(Dims)
 * to get the proper dimension number in the subexpression.
 */
template <unsigned NDim, unsigned Dim, typename T, typename... Dims>
detail::enable_if_t<(Dim >= NDim-sizeof...(Dims))>
increment(array_expr<T, Dims...>& expr)
{
    increment(expr, std::get<Dim-(NDim-sizeof...(Dims))>(expr.dims));
}

template <unsigned NDim, unsigned Dim, typename T, typename... Dims>
detail::enable_if_t<(Dim >= NDim-sizeof...(Dims))>
decrement(array_expr<T, Dims...>& expr)
{
    decrement(expr, std::get<Dim-(NDim-sizeof...(Dims))>(expr.dims));
}

/*
 * Lastly, for scalars do nothing since they are implicitly broadcast (this
 * only happens when assigning directly to a scalar as any other scalars are
 * absorbed into expression nodes).
 */
template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<std::is_arithmetic<Expr>::value>
increment(const Expr& expr) {}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<std::is_arithmetic<Expr>::value>
decrement(const Expr& expr) {}

/*
 * For binary and unary subexpressions, increment/decrement their children.
 */
template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<is_binary_expression<Expr>::value>
increment(Expr& expr)
{
    increment<NDim, Dim>(expr.first);
    increment<NDim, Dim>(expr.second);
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<is_binary_expression<Expr>::value>
decrement(Expr& expr)
{
    decrement<NDim, Dim>(expr.first);
    decrement<NDim, Dim>(expr.second);
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<is_unary_expression<Expr>::value>
increment(Expr& expr)
{
    increment<NDim, Dim>(expr.expr);
}

template <unsigned NDim, unsigned Dim, typename Expr>
detail::enable_if_t<is_unary_expression<Expr>::value>
decrement(Expr& expr)
{
    decrement<NDim, Dim>(expr.expr);
}

template <unsigned NDim, unsigned Dim=1>
struct assign_expr_loop;

template <unsigned NDim>
struct assign_expr_loop<NDim, NDim>
{
    template <typename LHS, typename RHS>
    void operator()(LHS& lhs, RHS& rhs, const std::array<idx_type, NDim>& len) const
    {
        for (idx_type i = 0;i < len[NDim-1];i++)
        {
            eval(lhs) = eval(rhs);

            increment<NDim, NDim-1>(lhs);
            increment<NDim, NDim-1>(rhs);
        }

        decrement<NDim, NDim-1>(lhs);
        decrement<NDim, NDim-1>(rhs);
    }
};

template <unsigned NDim, unsigned Dim>
struct assign_expr_loop
{
    template <typename LHS, typename RHS>
    void operator()(LHS& lhs, RHS& rhs, const std::array<idx_type, NDim>& len) const
    {
        assign_expr_loop<NDim, Dim+1> next_loop;

        for (idx_type i = 0;i < len[Dim-1];i++)
        {
            next_loop(lhs, rhs, len);

            increment<NDim, Dim-1>(lhs);
            increment<NDim, Dim-1>(rhs);
        }

        decrement<NDim, Dim-1>(lhs);
        decrement<NDim, Dim-1>(rhs);
    }
};

template <unsigned NDim, unsigned Dim=1>
struct assign_expr_loop_vec_row_major;

template <unsigned NDim>
struct assign_expr_loop_vec_row_major<NDim, NDim>
{
    template <typename LHS, typename RHS>
    void operator()(LHS& lhs, RHS& rhs, const std::array<idx_type, NDim>& len) const
    {
        for (idx_type i = 0;i < len[NDim-1];i++)
        {
            eval_at<NDim, NDim-1>(lhs, i) = eval_at<NDim, NDim-1>(rhs, i);
        }
    }
};

template <unsigned NDim, unsigned Dim>
struct assign_expr_loop_vec_row_major
{
    template <typename LHS, typename RHS>
    void operator()(LHS& lhs, RHS& rhs, const std::array<idx_type, NDim>& len) const
    {
        assign_expr_loop_vec_row_major<NDim, Dim+1> next_loop;

        for (idx_type i = 0;i < len[Dim-1];i++)
        {
            next_loop(lhs, rhs, len);

            increment<NDim, Dim-1>(lhs);
            increment<NDim, Dim-1>(rhs);
        }

        decrement<NDim, Dim-1>(lhs);
        decrement<NDim, Dim-1>(rhs);
    }
};

template <unsigned NDim, unsigned Dim=NDim-1>
struct assign_expr_loop_vec_col_major;

template <unsigned NDim>
struct assign_expr_loop_vec_col_major<NDim, 0>
{
    template <typename LHS, typename RHS>
    void operator()(LHS& lhs, RHS& rhs, const std::array<idx_type, NDim>& len) const
    {
        for (idx_type i = 0;i < len[0];i++)
        {
            eval_at<NDim, 0>(lhs, i) = eval_at<NDim, 0>(rhs, i);
        }
    }
};

template <unsigned NDim, unsigned Dim>
struct assign_expr_loop_vec_col_major
{
    template <typename LHS, typename RHS>
    void operator()(LHS& lhs, RHS& rhs, const std::array<idx_type, NDim>& len) const
    {
        assign_expr_loop_vec_col_major<NDim, Dim-1> next_loop;

        for (idx_type i = 0;i < len[Dim];i++)
        {
            next_loop(lhs, rhs, len);

            increment<NDim, Dim>(lhs);
            increment<NDim, Dim>(rhs);
        }

        decrement<NDim, Dim>(lhs);
        decrement<NDim, Dim>(rhs);
    }
};

template <typename Array, typename Expr>
detail::enable_if_t<(is_array_expression<detail::decay_t<Array>>::value ||
                     is_marray<detail::decay_t<Array>>::value) &&
                    is_expression_arg_or_scalar<detail::decay_t<Expr>>::value>
assign_expr(Array&& array_, Expr&& expr_)
{
    typedef typename expression_type<detail::decay_t<Array>>::type array_type;
    typedef typename expression_type<detail::decay_t<Expr>>::type expr_type;

    static_assert(expr_dimension<array_type>::value >=
                  expr_dimension<expr_type>::value,
                  "Dimensionality of the expression must not exceed that of the target");

    constexpr unsigned NDim = expr_dimension<array_type>::value;

    auto array = make_expression(std::forward<Array>(array_));
    auto expr = make_expression(std::forward<Expr>(expr_));

    auto len = get_array_lengths(array);
    MARRAY_ASSERT(check_expr_lengths(expr, len));

    if (is_vectorizable<NDim, NDim-1>(array) &&
        is_vectorizable<NDim, NDim-1>(expr))
    {
        assign_expr_loop_vec_row_major<NDim>()(array, expr, len);
    }
    else if (is_vectorizable<NDim, 0>(array) &&
             is_vectorizable<NDim, 0>(expr))
    {
        assign_expr_loop_vec_col_major<NDim>()(array, expr, len);
    }
    else
    {
        assign_expr_loop<NDim>()(array, expr, len);
    }
}

}

#endif
