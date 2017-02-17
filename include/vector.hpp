#ifndef _MARRAY_VECTOR_HPP_
#define _MARRAY_VECTOR_HPP_

#include <complex>
#include "utility.hpp"

namespace MArray
{

template <typename T, typename=void>
struct vector_traits
{
    constexpr static unsigned vector_width = 1;
};

}

#if defined(__AVX512F__)

#include "vector_avx512.hpp"

#elif defined(__AVX__)

#include "vector_avx.hpp"

#elif defined(__SSE4_1__)

#include "vector_sse41.hpp"

#endif

#endif
