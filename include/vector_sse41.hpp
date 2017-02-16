#ifndef _MARRAY_VECTOR_SSE3_HPP_
#define _MARRAY_VECTOR_SSE3_HPP_

#include <x86intrin.h>
#include "vector.hpp"

namespace MArray
{

template <>
struct vector_traits<float>
{
    constexpr static unsigned vector_width = 4;
    constexpr static size_t alignment = 16;

    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128 v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128 v)
    {
        return _mm_cvtps_pd(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128 v)
    {
        return _mm_unpacklo_ps(v, _mm_setzero_ps());
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128 v)
    {
        __m128i i32 = _mm_cvtps_epi32(v);
        __m128i i16 = _mm_packs_epi32(i32, i32);
        return _mm_packs_epi16(i16, i16);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128 v)
    {
        __m128i i32 = _mm_cvtps_epi32(v);
        return _mm_packs_epi32(i32, i32);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128 v)
    {
        return _mm_cvtps_epi32(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128 v)
    {
        return _mm_set_epi64x((int64_t)v[0], (int64_t)v[1]);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned, __m128>
    load(const float* ptr)
    {
        return _mm_loadu_ps(ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned, __m128>
    load(const float* ptr)
    {
        return _mm_load_ps(ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2, __m128>
    load(const float* ptr)
    {
        return _mm_castpd_ps(_mm_load1_pd((double*)ptr));
    }

    __m128 load1(const float* ptr)
    {
        return _mm_load1_ps(ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned>
    store(__m128 v, float* ptr)
    {
        _mm_storeu_ps(ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned>
    store(__m128 v, float* ptr)
    {
        _mm_store_ps(ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2>
    store(__m128 v, float* ptr)
    {
        _mm_store_sd((double*)ptr, _mm_castpd_ps(v));
    }

    __m128 add(__m128 a, __m128 b)
    {
        return _mm_add_ps(a, b);
    }

    __m128 sub(__m128 a, __m128 b)
    {
        return _mm_sub_ps(a, b);
    }

    __m128 mul(__m128 a, __m128 b)
    {
        return _mm_mul_ps(a, b);
    }

    __m128 div(__m128 a, __m128 b)
    {
        return _mm_div_ps(a, b);
    }

    __m128 pow(__m128 a, __m128 b)
    {
        return _mm_set_ps(std::pow((float)a[0], (float)b[0]),
                          std::pow((float)a[1], (float)b[1]),
                          std::pow((float)a[2], (float)b[2]),
                          std::pow((float)a[3], (float)b[3]));
    }

    __m128 negate(__m128 a)
    {
        return _mm_xor_ps(a, _mm_set1_ps(-0.0f));
    }

    __m128 exp(__m128 a)
    {
        return _mm_set_ps(std::exp((float)a[0]),
                          std::exp((float)a[1]),
                          std::exp((float)a[2]),
                          std::exp((float)a[3]));
    }

    __m128 sqrt(__m128 a)
    {
        return _mm_sqrt_ps(a);
    }
};

template <>
struct vector_traits<double>
{
    constexpr static unsigned vector_width = 2;
    constexpr static size_t alignment = 16;

    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128d v)
    {
        return _mm_cvtpd_ps(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128d v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128d v)
    {
        return _mm_unpacklo_ps(_mm_cvtpd_ps(v), _mm_setzero_ps());
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128d v)
    {
        __m128i i32 = _mm_cvtpd_epi32(v);
        __m128i i16 = _mm_packs_epi32(i32, i32);
        return _mm_packs_epi16(i16, i16);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128d v)
    {
        __m128i i32 = _mm_cvtpd_epi32(v);
        return _mm_packs_epi32(i32, i32);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128d v)
    {
        return _mm_cvtpd_epi32(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128d v)
    {
        return _mm_set_epi64x((int64_t)v[0], (int64_t)v[1]);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && !Aligned, __m128d>
    load(const double* ptr)
    {
        return _mm_loadu_pd(ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && Aligned, __m128d>
    load(const double* ptr)
    {
        return _mm_load_pd(ptr);
    }

    __m128d load1(const double* ptr)
    {
        return _mm_load1_pd(ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && !Aligned>
    store(__m128d v, double* ptr)
    {
        _mm_storeu_pd(ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && Aligned>
    store(__m128d v, double* ptr)
    {
        _mm_store_pd(ptr, v);
    }

    __m128d add(__m128d a, __m128d b)
    {
        return _mm_add_pd(a, b);
    }

    __m128d sub(__m128d a, __m128d b)
    {
        return _mm_sub_pd(a, b);
    }

    __m128d mul(__m128d a, __m128d b)
    {
        return _mm_mul_pd(a, b);
    }

    __m128d div(__m128d a, __m128d b)
    {
        return _mm_div_pd(a, b);
    }

    __m128d pow(__m128d a, __m128d b)
    {
        return _mm_set_pd(std::pow((double)a[0], (double)b[0]),
                          std::pow((double)a[1], (double)b[1]));
    }

    __m128d negate(__m128d a)
    {
        return _mm_sub_pd(_mm_set1_pd(0.0), a);
    }

    __m128d exp(__m128d a)
    {
        return _mm_set_pd(std::exp((double)a[0]),
                          std::exp((double)a[1]));
    }

    __m128d sqrt(__m128d a)
    {
        return _mm_sqrt_pd(a);
    }
};

template <>
struct vector_traits<std::complex<float>>
{
    constexpr static unsigned vector_width = 2;
    constexpr static size_t alignment = 16;

    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128 v)
    {
        return _mm_shuffle_ps(v, v, _MM_SHUFFLE(3,1,2,0));
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128 v)
    {
        return _mm_cvtps_pd(convert<float>(v));
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128 v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128 v)
    {
        __m128i i32 = _mm_cvtps_epi32(convert<float>(v));
        __m128i i16 = _mm_packs_epi32(i32, i32);
        return _mm_packs_epi16(i16, i16);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128 v)
    {
        __m128i i32 = _mm_cvtps_epi32(convert<float>(v));
        return _mm_packs_epi32(i32, i32);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128 v)
    {
        return _mm_cvtps_epi32(convert<float>(v));
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128 v)
    {
        return _mm_set_epi64x((int64_t)v[0], (int64_t)v[2]);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && !Aligned, __m128>
    load(const std::complex<float>* ptr)
    {
        return _mm_loadu_ps((float*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && Aligned, __m128>
    load(const std::complex<float>* ptr)
    {
        return _mm_load_ps((float*)ptr);
    }

    __m128 load1(const std::complex<float>* ptr)
    {
        return _mm_castpd_ps(_mm_load1_pd((double*)ptr));
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && !Aligned>
    store(__m128 v, std::complex<float>* ptr)
    {
        _mm_storeu_ps((float*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2 && Aligned>
    store(__m128 v, std::complex<float>* ptr)
    {
        _mm_store_ps((float*)ptr, v);
    }

    __m128 add(__m128 a, __m128 b)
    {
        return _mm_add_ps(a, b);
    }

    __m128 sub(__m128 a, __m128 b)
    {
        return _mm_sub_ps(a, b);
    }

    __m128 mul(__m128 a, __m128 b)
    {
        __m128 ashuf = _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,3,0,1));
        __m128 breal = _mm_moveldup_ps(b);
        __m128 bimag = _mm_movehdup_ps(b);
        __m128 tmp1 = _mm_mul_ps(    a, breal); // tmp1 = (ar0*br0, ai0*br0, ar1*br1, ai1*br1)
        __m128 tmp2 = _mm_mul_ps(ashuf, bimag); // tmp2 = (ai0*bi0, ar0*bi0, ai1*bi1, ar1*bi1)
        return _mm_addsub_ps(tmp1, tmp2);
    }

    __m128 div(__m128 a, __m128 b)
    {
        __m128 bsqr = _mm_mul_ps(b, b);
        bsqr = _mm_hadd_ps(bsqr, bsqr);
        bsqr = _mm_moveldup_ps(bsqr); // bsqr = (|b0|^2, |b0|^2, |b1|^2, |b1|^2)

        __m128 ashuf = _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,3,0,1));
        __m128 breal = _mm_moveldup_ps(b);
        __m128 bimag = _mm_movehdup_ps(b);
        __m128 tmp1 = _mm_mul_ps(    a, breal); // tmp1 = ( ar0*br0,  ai0*br0,  ar1*br1,  ai1*br1)
        __m128 tmp2 = _mm_mul_ps(ashuf, bimag);
        tmp2 = _mm_xor_ps(tmp2, _mm_set1_ps(-0.0f)); // tmp2 = (-ai0*bi0, -ar0*bi0, -ai1*bi1, -ar1*bi1)
        __m128 abconj = _mm_addsub_ps(tmp1, tmp2);

        return _mm_div_ps(abconj, bsqr);
    }

    __m128 pow(__m128 a, __m128 b)
    {
        std::complex<float> a0((float)a[0], (float)a[1]);
        std::complex<float> a1((float)a[2], (float)a[3]);
        std::complex<float> b0((float)b[0], (float)b[1]);
        std::complex<float> b1((float)b[2], (float)b[3]);
        std::complex<float> c0 = std::pow(a0, b0);
        std::complex<float> c1 = std::pow(a1, b1);
        return _mm_set_ps(c0.real(), c0.imag(),
                          c1.real(), c1.imag());
    }

    __m128 negate(__m128 a)
    {
        return _mm_xor_ps(a, _mm_set1_ps(-0.0f));
    }

    __m128 exp(__m128 a)
    {
        std::complex<float> a0((float)a[0], (float)a[1]);
        std::complex<float> a1((float)a[2], (float)a[3]);
        std::complex<float> b0 = std::exp(a0);
        std::complex<float> b1 = std::exp(a1);
        return _mm_set_ps(b0.real(), b0.imag(),
                          b1.real(), b1.imag());
    }

    __m128 sqrt(__m128 a)
    {
        std::complex<float> a0((float)a[0], (float)a[1]);
        std::complex<float> a1((float)a[2], (float)a[3]);
        std::complex<float> b0 = std::sqrt(a0);
        std::complex<float> b1 = std::sqrt(a1);
        return _mm_set_ps(b0.real(), b0.imag(),
                          b1.real(), b1.imag());
    }
};

template <>
struct vector_traits<int8_t>
{
    constexpr static unsigned vector_width = 16;
    constexpr static size_t alignment = 16;

    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128i v)
    {
        return _mm_cvtepi32_ps(_mm_cvtepi8_epi32(v));
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128i v)
    {
        return _mm_cvtepi32_pd(_mm_cvtepi8_epi32(v));
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128i v)
    {
        return _mm_unpacklo_ps(_mm_cvtepi32_ps(_mm_cvtepi8_epi32(v)), _mm_setzero_ps());
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128i v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi8_epi16(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi8_epi32(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi8_epi64(v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && !Aligned, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_loadu_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && Aligned, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_load_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 8, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_set1_epi64x(*(int64_t*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_set1_epi32(*(int32_t*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_set1_epi16(*(int16_t*)ptr);
    }

    __m128i load1(const int8_t* ptr)
    {
        return _mm_set1_epi8(*ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && !Aligned>
    store(__m128i v, int8_t* ptr)
    {
        _mm_storeu_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && Aligned>
    store(__m128i v, int8_t* ptr)
    {
        _mm_store_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 8>
    store(__m128i v, int8_t* ptr)
    {
        _mm_storel_epi64((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4>
    store(__m128i v, int8_t* ptr)
    {
        *(int32_t*)ptr = _mm_extract_epi32(v, 0);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2>
    store(__m128i v, int8_t* ptr)
    {
        *(int16_t*)ptr = _mm_extract_epi16(v, 0);
    }

    __m128i add(__m128i a, __m128i b)
    {
        return _mm_add_epi8(a, b);
    }

    __m128i sub(__m128i a, __m128i b)
    {
        return _mm_sub_epi8(a, b);
    }

    __m128i mul(__m128i a, __m128i b)
    {
        __m128i lo = _mm_and_si128(_mm_mullo_epi16(a, b), _mm_set1_epi16(0xff));
        __m128i hi = _mm_mullo_epi16(_mm_srli_epi16(a, 8),_mm_srli_epi16(b, 8));
        return _mm_or_si128(_mm_slli_epi16(hi, 8), lo);
    }

    __m128i div(__m128i a, __m128i b)
    {
        return _mm_set_epi8((int8_t)_mm_extract_epi8(a, 0) /
                            (int8_t)_mm_extract_epi8(b, 0),
                            (int8_t)_mm_extract_epi8(a, 1) /
                            (int8_t)_mm_extract_epi8(b, 1),
                            (int8_t)_mm_extract_epi8(a, 2) /
                            (int8_t)_mm_extract_epi8(b, 2),
                            (int8_t)_mm_extract_epi8(a, 3) /
                            (int8_t)_mm_extract_epi8(b, 3),
                            (int8_t)_mm_extract_epi8(a, 4) /
                            (int8_t)_mm_extract_epi8(b, 4),
                            (int8_t)_mm_extract_epi8(a, 5) /
                            (int8_t)_mm_extract_epi8(b, 5),
                            (int8_t)_mm_extract_epi8(a, 6) /
                            (int8_t)_mm_extract_epi8(b, 6),
                            (int8_t)_mm_extract_epi8(a, 7) /
                            (int8_t)_mm_extract_epi8(b, 7),
                            (int8_t)_mm_extract_epi8(a, 8) /
                            (int8_t)_mm_extract_epi8(b, 8),
                            (int8_t)_mm_extract_epi8(a, 9) /
                            (int8_t)_mm_extract_epi8(b, 9),
                            (int8_t)_mm_extract_epi8(a,10) /
                            (int8_t)_mm_extract_epi8(b,10),
                            (int8_t)_mm_extract_epi8(a,11) /
                            (int8_t)_mm_extract_epi8(b,11),
                            (int8_t)_mm_extract_epi8(a,12) /
                            (int8_t)_mm_extract_epi8(b,12),
                            (int8_t)_mm_extract_epi8(a,13) /
                            (int8_t)_mm_extract_epi8(b,13),
                            (int8_t)_mm_extract_epi8(a,14) /
                            (int8_t)_mm_extract_epi8(b,14),
                            (int8_t)_mm_extract_epi8(a,15) /
                            (int8_t)_mm_extract_epi8(b,15));
    }

    __m128i pow(__m128i a, __m128i b)
    {
        return _mm_set_epi8((int8_t)std::pow(_mm_extract_epi8(a, 0),
                                             _mm_extract_epi8(b, 0)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 1),
                                             _mm_extract_epi8(b, 1)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 2),
                                             _mm_extract_epi8(b, 2)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 3),
                                             _mm_extract_epi8(b, 3)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 4),
                                             _mm_extract_epi8(b, 4)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 5),
                                             _mm_extract_epi8(b, 5)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 6),
                                             _mm_extract_epi8(b, 6)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 7),
                                             _mm_extract_epi8(b, 7)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 8),
                                             _mm_extract_epi8(b, 8)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 9),
                                             _mm_extract_epi8(b, 9)),
                            (int8_t)std::pow(_mm_extract_epi8(a,10),
                                             _mm_extract_epi8(b,10)),
                            (int8_t)std::pow(_mm_extract_epi8(a,11),
                                             _mm_extract_epi8(b,11)),
                            (int8_t)std::pow(_mm_extract_epi8(a,12),
                                             _mm_extract_epi8(b,12)),
                            (int8_t)std::pow(_mm_extract_epi8(a,13),
                                             _mm_extract_epi8(b,13)),
                            (int8_t)std::pow(_mm_extract_epi8(a,14),
                                             _mm_extract_epi8(b,14)),
                            (int8_t)std::pow(_mm_extract_epi8(a,15),
                                             _mm_extract_epi8(b,15)));                                                         _mm_extract_epi8(b, 3)));
    }

    __m128i negate(__m128i a)
    {
        return _mm_sign_epi8(a, _mm_set1_epi8(0x80));
    }

    __m128i exp(__m128i a)
    {
        return _mm_set_epi8((int8_t)std::exp(_mm_extract_epi8(a, 0)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 1)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 2)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 3)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 4)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 5)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 6)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 7)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 8)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 9)),
                            (int8_t)std::exp(_mm_extract_epi8(a,10)),
                            (int8_t)std::exp(_mm_extract_epi8(a,11)),
                            (int8_t)std::exp(_mm_extract_epi8(a,12)),
                            (int8_t)std::exp(_mm_extract_epi8(a,13)),
                            (int8_t)std::exp(_mm_extract_epi8(a,14)),
                            (int8_t)std::exp(_mm_extract_epi8(a,15)));
    }

    __m128i sqrt(__m128i a)
    {
        return _mm_set_epi8((int8_t)std::sqrt(_mm_extract_epi8(a, 0)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 1)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 2)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 3)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 4)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 5)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 6)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 7)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 8)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 9)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,10)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,11)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,12)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,13)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,14)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,15)));
    }
};

template <>
struct vector_traits<uint8_t> : vector_traits<int8_t>
{
    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128i v)
    {
        return _mm_cvtepi32_ps(_mm_cvtepi8_epi32(v));
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128i v)
    {
        return _mm_cvtepi32_pd(_mm_cvtepi8_epi32(v));
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128i v)
    {
        return _mm_unpacklo_ps(_mm_cvtepi32_ps(_mm_cvtepi8_epi32(v)), _mm_setzero_ps());
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128i v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi8_epi16(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi8_epi32(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi8_epi64(v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && !Aligned, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_loadu_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && Aligned, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_load_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 8, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_set1_epi64x(*(int64_t*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_set1_epi32(*(int32_t*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2, __m128i>
    load(const int8_t* ptr)
    {
        return _mm_set1_epi16(*(int16_t*)ptr);
    }

    __m128i load1(const int8_t* ptr)
    {
        return _mm_set1_epi8(*ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && !Aligned>
    store(__m128i v, int8_t* ptr)
    {
        _mm_storeu_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 16 && Aligned>
    store(__m128i v, int8_t* ptr)
    {
        _mm_store_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 8>
    store(__m128i v, int8_t* ptr)
    {
        _mm_storel_epi64((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4>
    store(__m128i v, int8_t* ptr)
    {
        *(int32_t*)ptr = _mm_extract_epi32(v, 0);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2>
    store(__m128i v, int8_t* ptr)
    {
        *(int16_t*)ptr = _mm_extract_epi16(v, 0);
    }

    __m128i add(__m128i a, __m128i b)
    {
        return _mm_add_epi8(a, b);
    }

    __m128i sub(__m128i a, __m128i b)
    {
        return _mm_sub_epi8(a, b);
    }

    __m128i mul(__m128i a, __m128i b)
    {
        //TODO
    }

    __m128i div(__m128i a, __m128i b)
    {
        return _mm_set_epi8((int8_t)_mm_extract_epi8(a, 0) /
                            (int8_t)_mm_extract_epi8(b, 0),
                            (int8_t)_mm_extract_epi8(a, 1) /
                            (int8_t)_mm_extract_epi8(b, 1),
                            (int8_t)_mm_extract_epi8(a, 2) /
                            (int8_t)_mm_extract_epi8(b, 2),
                            (int8_t)_mm_extract_epi8(a, 3) /
                            (int8_t)_mm_extract_epi8(b, 3),
                            (int8_t)_mm_extract_epi8(a, 4) /
                            (int8_t)_mm_extract_epi8(b, 4),
                            (int8_t)_mm_extract_epi8(a, 5) /
                            (int8_t)_mm_extract_epi8(b, 5),
                            (int8_t)_mm_extract_epi8(a, 6) /
                            (int8_t)_mm_extract_epi8(b, 6),
                            (int8_t)_mm_extract_epi8(a, 7) /
                            (int8_t)_mm_extract_epi8(b, 7),
                            (int8_t)_mm_extract_epi8(a, 8) /
                            (int8_t)_mm_extract_epi8(b, 8),
                            (int8_t)_mm_extract_epi8(a, 9) /
                            (int8_t)_mm_extract_epi8(b, 9),
                            (int8_t)_mm_extract_epi8(a,10) /
                            (int8_t)_mm_extract_epi8(b,10),
                            (int8_t)_mm_extract_epi8(a,11) /
                            (int8_t)_mm_extract_epi8(b,11),
                            (int8_t)_mm_extract_epi8(a,12) /
                            (int8_t)_mm_extract_epi8(b,12),
                            (int8_t)_mm_extract_epi8(a,13) /
                            (int8_t)_mm_extract_epi8(b,13),
                            (int8_t)_mm_extract_epi8(a,14) /
                            (int8_t)_mm_extract_epi8(b,14),
                            (int8_t)_mm_extract_epi8(a,15) /
                            (int8_t)_mm_extract_epi8(b,15));
    }

    __m128i pow(__m128i a, __m128i b)
    {
        return _mm_set_epi8((int8_t)std::pow(_mm_extract_epi8(a, 0),
                                             _mm_extract_epi8(b, 0)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 1),
                                             _mm_extract_epi8(b, 1)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 2),
                                             _mm_extract_epi8(b, 2)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 3),
                                             _mm_extract_epi8(b, 3)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 4),
                                             _mm_extract_epi8(b, 4)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 5),
                                             _mm_extract_epi8(b, 5)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 6),
                                             _mm_extract_epi8(b, 6)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 7),
                                             _mm_extract_epi8(b, 7)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 8),
                                             _mm_extract_epi8(b, 8)),
                            (int8_t)std::pow(_mm_extract_epi8(a, 9),
                                             _mm_extract_epi8(b, 9)),
                            (int8_t)std::pow(_mm_extract_epi8(a,10),
                                             _mm_extract_epi8(b,10)),
                            (int8_t)std::pow(_mm_extract_epi8(a,11),
                                             _mm_extract_epi8(b,11)),
                            (int8_t)std::pow(_mm_extract_epi8(a,12),
                                             _mm_extract_epi8(b,12)),
                            (int8_t)std::pow(_mm_extract_epi8(a,13),
                                             _mm_extract_epi8(b,13)),
                            (int8_t)std::pow(_mm_extract_epi8(a,14),
                                             _mm_extract_epi8(b,14)),
                            (int8_t)std::pow(_mm_extract_epi8(a,15),
                                             _mm_extract_epi8(b,15)));                                                         _mm_extract_epi8(b, 3)));
    }

    __m128i negate(__m128i a)
    {
        return _mm_sign_epi8(a, _mm_set1_epi8(0x80));
    }

    __m128i exp(__m128i a)
    {
        return _mm_set_epi8((int8_t)std::exp(_mm_extract_epi8(a, 0)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 1)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 2)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 3)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 4)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 5)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 6)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 7)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 8)),
                            (int8_t)std::exp(_mm_extract_epi8(a, 9)),
                            (int8_t)std::exp(_mm_extract_epi8(a,10)),
                            (int8_t)std::exp(_mm_extract_epi8(a,11)),
                            (int8_t)std::exp(_mm_extract_epi8(a,12)),
                            (int8_t)std::exp(_mm_extract_epi8(a,13)),
                            (int8_t)std::exp(_mm_extract_epi8(a,14)),
                            (int8_t)std::exp(_mm_extract_epi8(a,15)));
    }

    __m128i sqrt(__m128i a)
    {
        return _mm_set_epi8((int8_t)std::sqrt(_mm_extract_epi8(a, 0)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 1)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 2)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 3)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 4)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 5)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 6)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 7)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 8)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a, 9)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,10)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,11)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,12)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,13)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,14)),
                            (int8_t)std::sqrt(_mm_extract_epi8(a,15)));
    }
};

template <>
struct vector_traits<int32_t>
{
    constexpr static unsigned vector_width = 4;
    constexpr static size_t alignment = 16;

    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128i v)
    {
        return _mm_cvtepi32_ps(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128i v)
    {
        return _mm_cvtepi32_pd(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128i v)
    {
        return _mm_unpacklo_ps(_mm_cvtepi32_ps(v), _mm_setzero_ps());
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128i v)
    {
        __m128i i16 = _mm_packs_epi32(v, v);
        return _mm_packs_epi16(i16, i16);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_packs_epi32(v, v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128i v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi32_epi64(v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_loadu_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_load_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_set1_epi64x(*(int64_t*)ptr);
    }

    __m128i load1(const int32_t* ptr)
    {
        return _mm_set1_epi32(*ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned>
    store(__m128i v, int32_t* ptr)
    {
        _mm_storeu_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned>
    store(__m128i v, int32_t* ptr)
    {
        _mm_store_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2>
    store(__m128i v, int32_t* ptr)
    {
        _mm_storel_epi64((__m128i*)ptr, v);
    }

    __m128i add(__m128i a, __m128i b)
    {
        return _mm_add_epi32(a, b);
    }

    __m128i sub(__m128i a, __m128i b)
    {
        return _mm_sub_epi32(a, b);
    }

    __m128i mul(__m128i a, __m128i b)
    {
        return _mm_mullo_epi32(a, b);
    }

    __m128i div(__m128i a, __m128i b)
    {
        return _mm_set_epi32((int32_t)_mm_extract_epi32(a, 0) /
                             (int32_t)_mm_extract_epi32(b, 0),
                             (int32_t)_mm_extract_epi32(a, 1) /
                             (int32_t)_mm_extract_epi32(b, 1),
                             (int32_t)_mm_extract_epi32(a, 2) /
                             (int32_t)_mm_extract_epi32(b, 2),
                             (int32_t)_mm_extract_epi32(a, 3) /
                             (int32_t)_mm_extract_epi32(b, 3));
    }

    __m128i pow(__m128i a, __m128i b)
    {
        return _mm_set_epi32((int32_t)std::pow(_mm_extract_epi32(a, 0),
                                               _mm_extract_epi32(b, 0)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 1),
                                               _mm_extract_epi32(b, 1)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 2),
                                               _mm_extract_epi32(b, 2)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 3),
                                               _mm_extract_epi32(b, 3)));
    }

    __m128i negate(__m128i a)
    {
        return _mm_sign_epi32(a, _mm_set1_epi32(0x80000000));
    }

    __m128i exp(__m128i a)
    {
        return _mm_set_epi32((int32_t)std::exp(_mm_extract_epi32(a, 0)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 1)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 2)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 3)));
    }

    __m128i sqrt(__m128i a)
    {
        return _mm_set_epi32((int32_t)std::sqrt(_mm_extract_epi32(a, 0)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 1)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 2)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 3)));
    }
};

template <>
struct vector_traits<int32_t>
{
    constexpr static unsigned vector_width = 4;
    constexpr static size_t alignment = 16;

    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128i v)
    {
        return _mm_cvtepi32_ps(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128i v)
    {
        return _mm_cvtepi32_pd(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128i v)
    {
        return _mm_unpacklo_ps(_mm_cvtepi32_ps(v), _mm_setzero_ps());
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128i v)
    {
        __m128i i16 = _mm_packs_epi32(v, v);
        return _mm_packs_epi16(i16, i16);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_packs_epi32(v, v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128i v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi32_epi64(v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_loadu_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_load_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_set1_epi64x(*(int64_t*)ptr);
    }

    __m128i load1(const int32_t* ptr)
    {
        return _mm_set1_epi32(*ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned>
    store(__m128i v, int32_t* ptr)
    {
        _mm_storeu_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned>
    store(__m128i v, int32_t* ptr)
    {
        _mm_store_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2>
    store(__m128i v, int32_t* ptr)
    {
        _mm_storel_epi64((__m128i*)ptr, v);
    }

    __m128i add(__m128i a, __m128i b)
    {
        return _mm_add_epi32(a, b);
    }

    __m128i sub(__m128i a, __m128i b)
    {
        return _mm_sub_epi32(a, b);
    }

    __m128i mul(__m128i a, __m128i b)
    {
        return _mm_mullo_epi32(a, b);
    }

    __m128i div(__m128i a, __m128i b)
    {
        return _mm_set_epi32((int32_t)_mm_extract_epi32(a, 0) /
                             (int32_t)_mm_extract_epi32(b, 0),
                             (int32_t)_mm_extract_epi32(a, 1) /
                             (int32_t)_mm_extract_epi32(b, 1),
                             (int32_t)_mm_extract_epi32(a, 2) /
                             (int32_t)_mm_extract_epi32(b, 2),
                             (int32_t)_mm_extract_epi32(a, 3) /
                             (int32_t)_mm_extract_epi32(b, 3));
    }

    __m128i pow(__m128i a, __m128i b)
    {
        return _mm_set_epi32((int32_t)std::pow(_mm_extract_epi32(a, 0),
                                               _mm_extract_epi32(b, 0)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 1),
                                               _mm_extract_epi32(b, 1)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 2),
                                               _mm_extract_epi32(b, 2)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 3),
                                               _mm_extract_epi32(b, 3)));
    }

    __m128i negate(__m128i a)
    {
        return _mm_sign_epi32(a, _mm_set1_epi32(0x80000000));
    }

    __m128i exp(__m128i a)
    {
        return _mm_set_epi32((int32_t)std::exp(_mm_extract_epi32(a, 0)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 1)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 2)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 3)));
    }

    __m128i sqrt(__m128i a)
    {
        return _mm_set_epi32((int32_t)std::sqrt(_mm_extract_epi32(a, 0)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 1)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 2)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 3)));
    }
};

template <>
struct vector_traits<int32_t>
{
    constexpr static unsigned vector_width = 4;
    constexpr static size_t alignment = 16;

    template <typename T>
    detail::enable_if_t<std::is_same<T,float>::value, __m128>
    convert(__m128i v)
    {
        return _mm_cvtepi32_ps(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,double>::value, __m128d>
    convert(__m128i v)
    {
        return _mm_cvtepi32_pd(v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,std::complex<float>>::value, __m128>
    convert(__m128i v)
    {
        return _mm_unpacklo_ps(_mm_cvtepi32_ps(v), _mm_setzero_ps());
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int8_t>::value ||
                        std::is_same<T,uint8_t>::value, __m128i>
    convert(__m128i v)
    {
        __m128i i16 = _mm_packs_epi32(v, v);
        return _mm_packs_epi16(i16, i16);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int16_t>::value ||
                        std::is_same<T,uint16_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_packs_epi32(v, v);
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int32_t>::value ||
                        std::is_same<T,uint32_t>::value, __m128i>
    convert(__m128i v)
    {
        return v;
    }

    template <typename T>
    detail::enable_if_t<std::is_same<T,int64_t>::value ||
                        std::is_same<T,uint64_t>::value, __m128i>
    convert(__m128i v)
    {
        return _mm_cvtepi32_epi64(v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_loadu_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_load_si128((__m128i*)ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2, __m128i>
    load(const int32_t* ptr)
    {
        return _mm_set1_epi64x(*(int64_t*)ptr);
    }

    __m128i load1(const int32_t* ptr)
    {
        return _mm_set1_epi32(*ptr);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && !Aligned>
    store(__m128i v, int32_t* ptr)
    {
        _mm_storeu_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 4 && Aligned>
    store(__m128i v, int32_t* ptr)
    {
        _mm_store_si128((__m128i*)ptr, v);
    }

    template <unsigned Width, bool Aligned>
    detail::enable_if_t<Width == 2>
    store(__m128i v, int32_t* ptr)
    {
        _mm_storel_epi64((__m128i*)ptr, v);
    }

    __m128i add(__m128i a, __m128i b)
    {
        return _mm_add_epi32(a, b);
    }

    __m128i sub(__m128i a, __m128i b)
    {
        return _mm_sub_epi32(a, b);
    }

    __m128i mul(__m128i a, __m128i b)
    {
        return _mm_mullo_epi32(a, b);
    }

    __m128i div(__m128i a, __m128i b)
    {
        return _mm_set_epi32((int32_t)_mm_extract_epi32(a, 0) /
                             (int32_t)_mm_extract_epi32(b, 0),
                             (int32_t)_mm_extract_epi32(a, 1) /
                             (int32_t)_mm_extract_epi32(b, 1),
                             (int32_t)_mm_extract_epi32(a, 2) /
                             (int32_t)_mm_extract_epi32(b, 2),
                             (int32_t)_mm_extract_epi32(a, 3) /
                             (int32_t)_mm_extract_epi32(b, 3));
    }

    __m128i pow(__m128i a, __m128i b)
    {
        return _mm_set_epi32((int32_t)std::pow(_mm_extract_epi32(a, 0),
                                               _mm_extract_epi32(b, 0)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 1),
                                               _mm_extract_epi32(b, 1)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 2),
                                               _mm_extract_epi32(b, 2)),
                             (int32_t)std::pow(_mm_extract_epi32(a, 3),
                                               _mm_extract_epi32(b, 3)));
    }

    __m128i negate(__m128i a)
    {
        return _mm_sign_epi32(a, _mm_set1_epi32(0x80000000));
    }

    __m128i exp(__m128i a)
    {
        return _mm_set_epi32((int32_t)std::exp(_mm_extract_epi32(a, 0)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 1)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 2)),
                             (int32_t)std::exp(_mm_extract_epi32(a, 3)));
    }

    __m128i sqrt(__m128i a)
    {
        return _mm_set_epi32((int32_t)std::sqrt(_mm_extract_epi32(a, 0)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 1)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 2)),
                             (int32_t)std::sqrt(_mm_extract_epi32(a, 3)));
    }
};

}

#endif
