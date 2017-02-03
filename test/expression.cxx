#include "gtest/gtest.h"
#include "expression.hpp"

using namespace std;
using namespace MArray;

TEST(expression, assign)
{

}

#if 0
struct bcast_dim {};

struct slice_dim
{
    idx_type len;
    stride_type stride;

    slice_dim(idx_type len, stride_type stride)
    : len(len), stride(stride) {}
};

template <typename T, typename... Dims>
struct array_expr
{
    typedef T& result_type;

    T* data;
    std::tuple<Dims...> dims;

    array_expr(T* data, const Dims&... dims)
    : data(data), dims{dims...} {}

    result_type eval() const
    {
        return *data;
    }

    result_type eval_at(idx_type i) const
    {
        return data[i];
    }
};

template <typename Expr, typename=void> struct is_array_expression;

template <typename Expr, typename=void> struct is_expression;

template <typename Expr, typename=void> struct is_unary_expression;

template <typename Expr, typename=void> struct is_binary_expression;

template <typename Expr, typename=void> struct expr_result_type;

template <typename Expr>
struct expr_result_type<Expr, detail::enable_if_t<std::is_arithmetic<Expr>::value>
{
    typedef detail::decay_t<Expr> type;
};

template <typename Expr>
struct expr_result_type<Expr, detail::enable_if_t<is_expression<Expr>::value>
{
    typedef typename detail::decay_t<Expr>::result_type type;
};

template <typename Expr>
detail::enable_if_t<is_expression<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval(Expr&& expr);
{
    return expr.eval();
}

template <typename Expr>
detail::enable_if_t<std::is_arithmetic<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval(Expr&& expr);
{
    return expr;
}

template <typename Expr>
detail::enable_if_t<is_expression<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval_at(Expr&& expr, idx_type i);
{
    return expr.eval_at(i);
}

template <typename Expr>
detail::enable_if_t<std::is_arithmetic<detail::decay_t<Expr>>::value,
                    typename expr_result_type<detail::decay_t<Expr>>::type>
eval_at(Expr&& expr, idx_type i);
{
    return expr;
}

template <typename LHS, typename RHS>
struct add_expr
{
    typedef LHS first_type;
    typedef RHS second_type;
    typedef decltype(std::declval<typename expr_result_type<LHS>::type>() +
                     std::declval<typename expr_result_type<RHS>::type>()) result_type;

    LHS first;
    RHS second;

    add_expr(const LHS& first, const RHS& second)
    : first(first), second(second) {}

    result_type eval() const
    {
        return eval(first) + eval(second);
    }

    result_type eval_at(idx_type i) const
    {
        return eval_at(first, i) + eval_at(second, i);
    }
};

template <typename LHS, typename RHS>
struct sub_expr
{
    typedef LHS first_type;
    typedef RHS second_type;
    typedef decltype(std::declval<typename expr_result_type<LHS>::type>() -
                     std::declval<typename expr_result_type<RHS>::type>()) result_type;

    LHS first;
    RHS second;

    sub_expr(const LHS& first, const RHS& second)
    : first(first), second(second) {}

    result_type eval() const
    {
        return eval(first) - eval(second);
    }

    result_type eval_at(idx_type i) const
    {
        return eval_at(first, i) - eval_at(second, i);
    }
};

template <typename LHS, typename RHS>
struct mul_expr
{
    typedef LHS first_type;
    typedef RHS second_type;
    typedef decltype(std::declval<typename expr_result_type<LHS>::type>() *
                     std::declval<typename expr_result_type<RHS>::type>()) result_type;

    LHS first;
    RHS second;

    mul_expr(const LHS& first, const RHS& second)
    : first(first), second(second) {}

    result_type eval() const
    {
        return eval(first) * eval(second);
    }

    result_type eval_at(idx_type i) const
    {
        return eval_at(first, i) * eval_at(second, i);
    }
};

template <typename LHS, typename RHS>
struct div_expr
{
    typedef LHS first_type;
    typedef RHS second_type;
    typedef decltype(std::declval<typename expr_result_type<LHS>::type>() /
                     std::declval<typename expr_result_type<RHS>::type>()) result_type;

    LHS first;
    RHS second;

    div_expr(const LHS& first, const RHS& second)
    : first(first), second(second) {}

    result_type eval() const
    {
        return eval(first) / eval(second);
    }

    result_type eval_at(idx_type i) const
    {
        return eval_at(first, i) / eval_at(second, i);
    }
};

template <typename Base, typename Exponent>
struct pow_expr
{
    typedef Base first_type;
    typedef Exponent second_type;
    typedef decltype(std::pow(std::declval<typename expr_result_type<Base>::type>(),
                              std::declval<typename expr_result_type<Exponent>::type>())) result_type;

    Base first;
    Exponent second;

    pow_expr(const Base& first, const Exponent& second)
    : first(first), second(second) {}

    result_type eval() const
    {
        return std::pow(eval(first), eval(second));
    }

    result_type eval_at(idx_type i) const
    {
        return std::pow(eval_at(first, i), eval_at(second, i));
    }
};

template <typename Expr>
struct negate_expr
{
    typedef Expr expr_type;
    typedef decltype(-std::declval<typename expr_result_type<Expr>::type>()) result_type;

    Expr expr;

    negate_expr(const Expr& expr) : expr(expr) {}

    result_type eval() const
    {
        return -eval(expr);
    }

    result_type eval_at(idx_type i) const
    {
        return -eval_at(expr, i);
    }
};

template <typename Expr>
struct exp_expr
{
    typedef Expr expr_type;
    typedef decltype(std::exp(std::declval<typename expr_result_type<Expr>::type>())) result_type;

    Expr expr;

    exp_expr(const Expr& expr) : expr(expr) {}

    result_type eval() const
    {
        return std::exp(eval(expr));
    }

    result_type eval_at(idx_type i) const
    {
        return std::exp(eval_at(expr, i));
    }
};

template <typename Expr>
struct sqrt_expr
{
    typedef Expr expr_type;
    typedef decltype(std::sqrt(std::declval<typename expr_result_type<Expr>::type>())) result_type;

    Expr expr;

    sqrt_expr(const Expr& expr) : expr(expr) {}

    result_type eval() const
    {
        return std::sqrt(eval(expr));
    }

    result_type eval_at(idx_type i) const
    {
        return std::sqrt(eval_at(expr, i));
    }
};

template <typename T>
struct is_array_expression : std::false_type {};

template <typename T, typename... Dims>
struct is_array_expression<array_expr<T, Dims...>> : std::true_type {};

template <typename T>
struct is_expression : std::false_type {};

template <typename T, typename... Dims>
struct is_expression<array_expr<T, Dims...>> : std::true_type {};

template <typename LHS, typename RHS>
struct is_expression<add_expr<LHS, RHS>> : std::true_type {};

template <typename LHS, typename RHS>
struct is_expression<sub_expr<LHS, RHS>> : std::true_type {};

template <typename LHS, typename RHS>
struct is_expression<mul_expr<LHS, RHS>> : std::true_type {};

template <typename LHS, typename RHS>
struct is_expression<div_expr<LHS, RHS>> : std::true_type {};

template <typename Base, typename Exponent>
struct is_expression<pow_expr<Base, Exponent>> : std::true_type {};

template <typename Expr>
struct is_expression<negate_expr<Expr>> : std::true_type {};

template <typename Expr>
struct is_expression<exp_expr<Expr>> : std::true_type {};

template <typename Expr>
struct is_expression<sqrt_expr<Expr>> : std::true_type {};

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
    detail::enable_if_t<is_expression<typename Expr::first_type>::value>>
    : std::true_type {};

template <typename Array, typename Dim>
struct array_expr_type_helper2;

template <typename T, typename... Dims, typename Dim>
struct array_expr_type_helper2<array_expr<T, Dims...>, Dim>
{
    typedef array_expr<T, Dims..., Dim> type;
};

template <typename T, unsigned NDim>
struct array_expr_type_helper;

template <typename T>
struct array_expr_type_helper<T, 1>
{
    typedef array_expr<T, slice_dim> type;
};

template <typename T, unsigned NDim>
struct array_expr_type_helper
{
    typedef typename array_expr_type_helper2<
        typename array_expr_type_helper<T, NDim-1>::type, slice_dim>::type type;
};

template <typename Expr, typename=void>
struct expression_type;

template <typename T, unsigned NDim, unsigned NIndexed, unsigned NSliced>
struct expression_type<marray_slice<T, NDim, NIndexed, NSliced>>
{
    typedef typename array_expr_type_helper<T, NDim-NIndexed>::type type;
};

template <typename T, unsigned NDim>
struct expression_type<marray_view<T, NDim>>
{
    typedef typename array_expr_type_helper<T, NDim>::type type;
};

template <typename T, unsigned NDim, typename Alloc>
struct expression_type<marray<T, NDim, Alloc>>
{
    typedef typename array_expr_type_helper<T, NDim>::type type;
};

template <typename Expr>
struct expression_type<Expr, detail::enable_if_t<is_expression<Expr>::value ||
                                                 std::is_arithmetic<Expr>::value>>
{
    typedef Expr type;
};

template <typename T, unsigned NDim, unsigned NIndexed, unsigned NSliced,
          size_t... I, size_t... J>
typename expression_type<marray_slice<T, NDim, NIndexed, NSliced>>
make_expression_helper(const marray_slice<T, NDim, NIndexed, NSliced>& x,
                       detail::integer_sequence<size_t, I...>,
                       detail::integer_sequence<size_t, J...>)
{
    return {x.data(),
            sliced_dim(x.slice_length<I>(),
                       x.slice_stride<I>())...,
            sliced_dim(x.base().length<NIndexed+NSliced+J>(),
                       x.base().stride<NIndexed+NSliced+J>())...};
}

template <typename T, unsigned NDim, size_t... I>
typename expression_type<marray_view<T, NDim>>
make_expression_helper(const marray_view<T, NDim>& x,
                       detail::integer_sequence<size_t, I...>)
{
    return {x.data(), sliced_dim(x.length<I>(), x.stride<I>())...};
}

template <typename T, unsigned NDim, typename Alloc, size_t... I>
typename expression_type<marray<T, NDim, Alloc>>
make_expression_helper(const marray<T, NDim, Alloc>& x,
                       detail::integer_sequence<size_t, I...>)
{
    return {x.data(), sliced_dim(x.length<I>(), x.stride<I>())...};
}

template <typename T, unsigned NDim, unsigned NIndexed, unsigned NSliced>
typename expression_type<marray_slice<T, NDim, NIndexed, NSliced>>
make_expression(const marray_slice<T, NDim, NIndexed, NSliced>& x)
{
    return make_expression_helper(x, detail::static_range<NSliced>(),
                                  detail::static_range<NDim-NIndexed-NSliced>());
}

template <typename T, unsigned NDim>
typename expression_type<marray_view<T, NDim>>
make_expression(const marray_view<T, NDim>& x)
{
    return make_expression_helper(x, detail::static_range<NDim>());
}

template <typename T, unsigned NDim, typename Alloc>
typename expression_type<marray<T, NDim, Alloc>>
make_expression(const marray<T, NDim, Alloc>& x)
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

template <typename T, unsigned NDim, unsigned NIndexed, unsigned NSliced>
struct is_marray<marray_slice<T, NDim, NIndexed, NSliced>> : std::true_type {};

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
                    add_expr<typename expression_type<LHS>::type,
                             typename expression_type<RHS>::type>>
operator+(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename LHS, typename RHS>
detail::enable_if_t<are_expression_args<LHS,RHS>::value,
                    sub_expr<typename expression_type<LHS>::type,
                             typename expression_type<RHS>::type>>
operator-(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename LHS, typename RHS>
detail::enable_if_t<are_expression_args<LHS,RHS>::value,
                    mul_expr<typename expression_type<LHS>::type,
                             typename expression_type<RHS>::type>>
operator*(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename LHS, typename RHS>
detail::enable_if_t<are_expression_args<LHS,RHS>::value,
                    div_expr<typename expression_type<LHS>::type,
                             typename expression_type<RHS>::type>>
operator/(const LHS& lhs, const RHS& rhs)
{
    return {make_expression(lhs), make_expression(rhs)};
}

template <typename Base, typename Exponent>
detail::enable_if_t<are_expression_args<Base,Exponent>::value,
                    pow_expr<typename expression_type<Base>::type,
                             typename expression_type<Exponent>::type>>
pow(const Base& base, const Exponent& exponent)
{
    return {make_expression(base), make_expression(exponent)};
}

template <typename Expr>
detail::enable_if_t<is_expression_arg<Expr>::value,
                    negate_expr<typename expression_type<Expr>::type>>
operator-(const Expr& expr)
{
    return {make_expression(expr)};
}

template <typename Expr>
detail::enable_if_t<is_expression_arg<Expr>::value,
                    exp_expr<typename expression_type<Expr>::type>>
exp(const Expr& expr)
{
    return {make_expression(expr)};
}

template <typename Expr>
detail::enable_if_t<is_expression_arg<Expr>::value,
                    sqrt_expr<typename expression_type<Expr>::type>>
sqrt(const Expr& expr)
{
    return {make_expression(expr)};
}

template <typename Expr, typename=void> struct expr_dimension;

template <typename T, typename... Dims>
struct expr_dimension<array_expr<T, Dims...>>
    : std::integral_constant<unsigned, sizeof...(Dims)> {};

template <typename Expr>
struct expr_dimension<Expr, detail::enable_if_t<std::is_arithmetic<Expr>::value>
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

idx_type get_array_length(const bcast_dim&)
{
    static_assert(false, "Broadcast dimensions cannot be written to");
    return 0;
}

idx_type get_array_length(const slice_dim& dim)
{
    return dim.len;
}

template <typename T, typename... Dims, size_t... I>
std::array<idx_type, sizeof...(I)>
get_array_lengths_helper(const array_expr<T, Dims...>& array,
                         detail::integer_sequence<size_t, I...>)
{
    return {get_array_length(std::get<I>(array.dims)...)};
}

template <typename T, typename... Dims>
std::array<idx_type, sizeof...(Dims)>
get_array_lengths(const array_expr<T, Dims...>& array)
{
    return get_array_lengths_helper(array, detail::static_range<sizeof...(Dims)>());
}

bool check_expr_length(const bcast_dim&, idx_type)
{
    return true;
}

bool check_expr_length(const slice_dim& dim, idx_type len)
{
    return len == dim.len;
}

template <typename T, typename... Dims, size_t NDim, size_t... I>
void check_expr_lengths_helper(const array_expr<T, Dims...>& array,
                               const std::array<idx_type, NDim>& len,
                               detail::integer_sequence<size_t, I...>)
{
    std::array<bool, sizeof...(Dims)> values =
        {check_expr_length(std::get<I>(array.dims),
                           len[NDim-sizeof...(Dims)+I])...};

    bool ret = true;
    for (value : values) ret = ret && value;
    return ret;
}

template <typename T, typename... Dims, size_t NDim>
void check_expr_lengths(const array_expr<T, Dims...>& array,
                        const std::array<idx_type, NDim>& len)
{
    check_expr_lengths_helper(array, len, detail::static_range<sizeof...(Dims)>());
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

template <typename Expr, size_t NDim>
detail::enable_if_t<std::is_arithmetic<Expr>::value,bool>
check_expr_lengths(const Expr& expr, const std::array<idx_type, NDim>& len)
{
    return true;
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
increment(array_expr<T, Dims...>& expr)
{
    increment(expr, std::get<Dim>(expr.dims));
}

template <unsigned NDim, unsigned Dim, typename T, typename... Dims>
detail::enable_if_t<(Dim < NDim-sizeof...(Dims))>
decrement(array_expr<T, Dims...>& expr)
{
    decrement(expr, std::get<Dim>(expr.dims));
}

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

template <unsigned NDim, unsigned Dim>
struct assign_expr_loop;

template <unsigned NDim>
struct assign_expr_loop<NDim, NDim>
{
    template <typename LHS, typename RHS>
    void operator()(LHS&& lhs, RHS&& rhs, const std::array<idx_type, NDim>& len) const
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
    void operator()(LHS&& lhs, RHS&& rhs, const std::array<idx_type, NDim>& len) const
    {
        assign_expr_loop<NDim, Dim+1> next_loop;

        for (idx_type i = 0;i < len[Dim-1];i++)
        {
            next_loop();

            increment<NDim, Dim-1>(lhs);
            increment<NDim, Dim-1>(rhs);
        }

        decrement<NDim, Dim-1>(lhs);
        decrement<NDim, Dim-1>(rhs);
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
    check_expr_lengths(expr, len);

    assign_expr_loop<NDim, 1>()(array, expr, len);
}

}

#endif
