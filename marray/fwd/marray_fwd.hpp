#ifndef MARRAY_MARRAY_FWD_HPP
#define MARRAY_MARRAY_FWD_HPP

#include <memory>

#include "../types.hpp"

namespace MArray
{
/**
 * Special value which indicates that the number of dimensions is not known at compile time.
 *
 * @ingroup constants
 */
#if MARRAY_DOXYGEN
constexpr int DYNAMIC;
#else
constexpr int DYNAMIC = -1;
#endif

/**
 * A partially-indexed tensor.
 *
 * This type cannot be constructed directly, but is returned by indexing a tensor
 * or tensor view. Provides a limited API for:
 *
 * - Further indexing, slicing, or broadcasting
 * - Obtaining the data pointer
 * - Creating a view from the currently-indexed portion
 * - Participating in element-wise expressions (including as the left-hand side)
 *
 * In the latter two cases, the resulting view or expression leaves any unindexed dimensions intact, i.e. it
 * is as if the remaining dimensions were indexed with `[slice::all]`.
 *
 * @ingroup classes
 */
template <typename Type, int NDim, int NIndexed, typename... Dims>
class marray_slice;

/**
 * Tensor base class.
 *
 * This class should not be used directly. Use @ref marray and @ref marray_view instead.
 *
 * @ingroup classes
 */
template <typename Type, int NDim, typename Derived, bool Owner>
class marray_base;

/**
 * A tensor (multi-dimensional array) view, which may either be mutable or immutable.
 *
 * @tparam Type     The type of the tensor elements. The view is immutable if this is const-qualified.
 *
 * @tparam NDim     The number of tensor dimensions, must be positive or DYNAMIC. Default is DYNAMIC.
 *
 * @ingroup classes
 */
template <typename Type, int NDim=DYNAMIC>
class marray_view;

/**
 * A tensor (multi-dimensional array) container.
 *
 * @tparam Type         The type of the tensor elements.
 *
 * @tparam NDim         The number of tensor dimensions, must be positive.
 *
 * @tparam Allocator    An allocator. If not specified, `std::allocator<Type>` is used.
 *
 * @ingroup classes
 */
template <typename Type, int NDim=DYNAMIC, typename Allocator=std::allocator<Type>>
class marray;

/**
 * Alias for a 1-dimensional tensor view.
 *
 * @tparam Type         The type of the tensor elements.
 *
 * @ingroup types
 */
template <typename Type> using row_view = marray_view<Type, 1>;

/**
 * Alias for a 1-dimensional tensor.
 *
 * @tparam Type         The type of the tensor elements.
 *
 * @tparam Allocator    An allocator. If not specified, `std::allocator<Type>` is used.
 *
 * @ingroup types
 */
template <typename Type, typename Allocator=std::allocator<Type>> using row = marray<Type, 1, Allocator>;

/**
 * Alias for a 2-dimensional tensor view.
 *
 * @tparam Type         The type of the tensor elements.
 *
 * @ingroup types
 */
template <typename Type> using matrix_view = marray_view<Type, 2>;

/**
 * Alias for a 2-dimensional tensor.
 *
 * @tparam Type         The type of the tensor elements.
 *
 * @tparam Allocator    An allocator. If not specified, `std::allocator<Type>` is used.
 *
 * @ingroup types
 */
template <typename Type, typename Allocator=std::allocator<Type>> using matrix = marray<Type, 2, Allocator>;

struct layout
{
    struct construct {};

    int type;

    constexpr explicit layout(int type, construct) : type(type) {}

    bool operator==(layout other) const { return type == other.type; }
    bool operator!=(layout other) const { return type != other.type; }
};

struct column_major_layout : layout { constexpr column_major_layout() : layout(0, construct{}) {} };

/**
 * Specifies that elements should be laid out in column-major order.
 *
 * @ingroup constants
 */
constexpr column_major_layout COLUMN_MAJOR;

struct row_major_layout : layout { constexpr row_major_layout() : layout(1, construct{}) {} };

/**
 * Specifies that elements should be laid out in row-major order.
 *
 * @ingroup constants
 */
constexpr row_major_layout ROW_MAJOR;

/**
 * The default layout for tensors (either row- or column-major).
 * This value is controlled by the macro @ref MARRAY_DEFAULT_LAYOUT.
 *
 * @ingroup constants
 */
constexpr decltype(MARRAY_DEFAULT_LAYOUT) DEFAULT_LAYOUT;

struct base
{
    struct construct {};

    int type;

    constexpr explicit base(int type, construct) : type(type) {}

    bool operator==(base other) const { return type == other.type; }
    bool operator!=(base other) const { return type != other.type; }
};

struct base_zero : base { constexpr base_zero() : base(0, construct{}) {} };

/**
 * Specifies that indices should start at 0.
 *
 * @ingroup constants
 */
constexpr base_zero BASE_ZERO;

struct base_one : base { constexpr base_one() : base(1, construct{}) {} };

/**
 * Specifies that indices should start at 1.
 *
 * @ingroup constants
 */
constexpr base_one BASE_ONE;

struct fortran_t
{
    operator layout() { return COLUMN_MAJOR; }
    operator base() { return BASE_ONE; }
};

/**
 * Specifies column-major layout and base-1 indices.
 *
 * @ingroup constants
 */
constexpr fortran_t FORTRAN;

/**
 * Specifies column-major layout and base-1 indices.
 *
 * @ingroup constants
 */
constexpr fortran_t MATLAB;

/**
 * The default base for indices (either 0 or 1). This value is controlled by the macro @ref MARRAY_DEFAULT_BASE.
 *
 * @ingroup constants
 */
constexpr decltype(MARRAY_DEFAULT_BASE) DEFAULT_BASE;

}

#endif //MARRAY_MARRAY_FWD_HPP