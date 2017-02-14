#ifndef _MARRAY_VECTOR_SSE3_HPP_
#define _MARRAY_VECTOR_SSE3_HPP_

#include <x86intrin.h>
#include "vector.hpp"

namespace MArray
{

    template <>
    struct vector_traits<float>
    {
        typedef float T;
        typedef __m128 vector_type;
        constexpr static unsigned vector_width = 4;
        constexpr static size_t alignment = 16;

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && !Aligned, vector_type>
        load(const T* ptr) { return _mm_loadu_ps(ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && Aligned, vector_type>
        load(const T* ptr) { return _mm_load_ps(ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,double>::value, vector_type>
        load(const T* ptr)
        {
            vector_type f0 = _mm_load_ss(ptr);
            vector_type f1 = _mm_load_ss(ptr+1);
            vector_type f01 = _mm_unpacklo_ps(f0, f1);
            return _mm_cvtps_pd(f01);
        }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,std::complex<float>>::value, vector_type>
        load(const T* ptr)
        {
            vector_type f0 = _mm_load_ss(ptr);
            vector_type f1 = _mm_load_ss(ptr+1);
            return _mm_movelh_ps(f0, f1);
        }

        template <typename cvt_type>
        detail::enable_if_t<std::is_same<cvt_type,T>::value, vector_type>
        load1(const T* ptr) { return _mm_load1_ps(ptr); }

        template <typename cvt_type>
        detail::enable_if_t<std::is_same<cvt_type,double>::value, vector_type>
        load1(const T* ptr)
        {
            vector_type f = _mm_load1_ps(ptr);
            return _mm_cvtps_pd(f);
        }

        template <typename cvt_type>
        detail::enable_if_t<std::is_same<cvt_type,std::complex<double>>::value, vector_type>
        load1(const T* ptr)
        {
            vector_type f0 = _mm_load_ss(ptr);
            return _mm_movelh_ps(f0, f0);
        }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && !Aligned>
        store(vector_type v, T* ptr) { _mm_storeu_ps(ptr, v); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && Aligned>
        store(vector_type v, T* ptr) { _mm_store_ps(ptr, v); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,double>::value>
        store(vector_type v, T* ptr)
        {
            vector_type f01 = _mm_cvtpd_ps(v);
            _mm_store_ss(ptr, f01);
            vector_type f1 = _mm_shuffle_ps(f01, f01, _MM_SHUFFLE(3,2,0,1));
            _mm_store_ss(ptr+1, f1);
        }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,std::complex<double>>::value>
        store(vector_type v, T* ptr)
        {
            _mm_store_ss(ptr, v);
            vector_type f1 = _mm_unpackhi_ps(v, v);
            _mm_store_ss(ptr+1, f1);
        }

        vector_type add(vector_type a, vector_type b) { return _mm_add_ps(a, b); }

        vector_type sub(vector_type a, vector_type b) { return _mm_sub_ps(a, b); }

        vector_type mul(vector_type a, vector_type b) { return _mm_mul_ps(a, b); }

        vector_type div(vector_type a, vector_type b) { return _mm_div_ps(a, b); }

        vector_type pow(vector_type a, vector_type b)
        {
            return _mm_set_ps(std::pow((float)a[0], (float)b[0]),
                              std::pow((float)a[1], (float)b[1]),
                              std::pow((float)a[2], (float)b[2]),
                              std::pow((float)a[3], (float)b[3]));
        }

        vector_type negate(vector_type a) { return _mm_sub_ps(_mm_set1_ps(0.0f), a); }

        vector_type exp(vector_type a)
        {
            return _mm_set_ps(std::exp((float)a[0]),
                              std::exp((float)a[1]),
                              std::exp((float)a[2]),
                              std::exp((float)a[3]));
        }

        vector_type sqrt(vector_type a) { return _mm_sqrt_ps(a); }
    };

    template <>
    struct vector_traits<double>
    {
        typedef double T;
        typedef __m128d vector_type;
        constexpr static unsigned vector_width = 2;
        constexpr static size_t alignment = 16;

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && !Aligned, vector_type>
        load(const T* ptr) { return _mm_loadu_pd(ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && Aligned, vector_type>
        load(const T* ptr) { return _mm_load_pd(ptr); }

        template <typename cvt_type>
        detail::enable_if_t<std::is_same<cvt_type,T>::value, vector_type>
        load1(const T* ptr) { return _mm_load1_pd(ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && !Aligned>
        store(vector_type v, T* ptr) { _mm_storeu_pd(ptr, v); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && Aligned>
        store(vector_type v, T* ptr) { _mm_store_pd(ptr, v); }

        vector_type add(vector_type a, vector_type b) { return _mm_add_pd(a, b); }

        vector_type sub(vector_type a, vector_type b) { return _mm_sub_pd(a, b); }

        vector_type mul(vector_type a, vector_type b) { return _mm_mul_pd(a, b); }

        vector_type div(vector_type a, vector_type b) { return _mm_div_pd(a, b); }

        vector_type pow(vector_type a, vector_type b)
        {
            return _mm_set_pd(std::pow((double)a[0], (double)b[0]),
                              std::pow((double)a[1], (double)b[1]));
        }

        vector_type negate(vector_type a) { return _mm_sub_pd(_mm_set1_pd(0.0), a); }

        vector_type exp(vector_type a)
        {
            return _mm_set_pd(std::exp((double)a[0]),
                              std::exp((double)a[1]));
        }

        vector_type sqrt(vector_type a) { return _mm_sqrt_pd(a); }
    };

    template <>
    struct vector_traits<std::complex<float>>
    {
        typedef std::complex<float> T;
        typedef __m128 vector_type;
        constexpr static unsigned vector_width = 2;
        constexpr static size_t alignment = 16;

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && !Aligned, vector_type>
        load(const T* ptr) { return _mm_loadu_ps((float*)ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && Aligned, vector_type>
        load(const T* ptr) { return _mm_load_ps((float*)ptr); }

        template <typename cvt_type>
        detail::enable_if_t<std::is_same<cvt_type,T>::value, vector_type>
        load1(const T* ptr) { return _mm_load1_ps((float*)ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && !Aligned>
        store(vector_type v, T* ptr) { _mm_storeu_ps((float*)ptr, v); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && Aligned>
        store(vector_type v, T* ptr) { _mm_store_ps((float*)ptr, v); }

        vector_type add(vector_type a, vector_type b) { return _mm_add_ps(a, b); }

        vector_type sub(vector_type a, vector_type b) { return _mm_sub_ps(a, b); }

        vector_type mul(vector_type a, vector_type b)
        {
            vector_type ashuf = _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,3,0,1));
            vector_type breal = _mm_moveldup_ps(b);
            vector_type bimag = _mm_movehdup_ps(b);
            vector_type tmp1 = _mm_mul_ps(    a, breal); // tmp1 = (ar0*br0, ai0*br0, ar1*br1, ai1*br1)
            vector_type tmp2 = _mm_mul_ps(ashuf, bimag); // tmp2 = (ai0*bi0, ar0*bi0, ai1*bi1, ar1*bi1)
            return _mm_addsub_ps(tmp1, tmp2);
        }

        vector_type div(vector_type a, vector_type b)
        {
            vector_type bsqr = _mm_mul_ps(b, b);
            bsqr = _mm_hadd_ps(bsqr, bsqr);
            bsqr = _mm_moveldup_ps(bsqr); // bsqr = (|b0|^2, |b0|^2, |b1|^2, |b1|^2)

            vector_type ashuf = _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,3,0,1));
            vector_type breal = _mm_moveldup_ps(b);
            vector_type bimag = _mm_movehdup_ps(b);
            vector_type tmp1 = _mm_mul_ps(    a, breal); // tmp1 = ( ar0*br0,  ai0*br0,  ar1*br1,  ai1*br1)
            vector_type tmp2 = _mm_mul_ps(ashuf, bimag);
            tmp2 = _mm_sub_ps(_mm_set1_ps(0.0f), tmp2); // tmp2 = (-ai0*bi0, -ar0*bi0, -ai1*bi1, -ar1*bi1)
            vector_type abconj = _mm_addsub_ps(tmp1, tmp2);

            return _mm_div_ps(abconj, bsqr);
        }

        vector_type pow(vector_type a, vector_type b)
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

        vector_type negate(vector_type a) { return _mm_sub_ps(_mm_set1_ps(0.0f), a); }

        vector_type exp(vector_type a)
        {
            std::complex<float> a0((float)a[0], (float)a[1]);
            std::complex<float> a1((float)a[2], (float)a[3]);
            std::complex<float> b0 = std::exp(a0);
            std::complex<float> b1 = std::exp(a1);
            return _mm_set_ps(b0.real(), b0.imag(),
                              b1.real(), b1.imag());
        }

        vector_type sqrt(vector_type a) { return _mm_sqrt_ps(a); }
    };

    template <>
    struct vector_traits<int32_t>
    {
        typedef int32_t T;
        typedef __m128i vector_type;
        constexpr static unsigned vector_width = 4;
        constexpr static size_t alignment = 16;

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<typename std::make_signed<cvt_type>::type,T>::value && !Aligned, vector_type>
        load(const T* ptr) { return _mm_loadu_si128((vector_type*)ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<typename std::make_signed<cvt_type>::type,T>::value && Aligned, vector_type>
        load(const T* ptr) { return _mm_load_si128((vector_type*)ptr); }

        template <typename cvt_type>
        detail::enable_if_t<std::is_same<typename std::make_signed<cvt_type>::type,T>::value, vector_type>
        load1(const T* ptr) { return _mm_set1_epi32(*ptr); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && !Aligned>
        store(vector_type v, T* ptr) { _mm_storeu_si128((vector_type*)ptr, v); }

        template <typename cvt_type, bool Aligned>
        detail::enable_if_t<std::is_same<cvt_type,T>::value && Aligned>
        store(vector_type v, T* ptr) { _mm_store_si128((vector_type*)ptr, v); }

        vector_type add(vector_type a, vector_type b) { return _mm_add_epi32(a, b); }

        vector_type sub(vector_type a, vector_type b) { return _mm_sub_epi32(a, b); }

        vector_type mul(vector_type a, vector_type b) { return _mm_mullo_epi32(a, b); }

        vector_type div(vector_type a, vector_type b)
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

        vector_type pow(vector_type a, vector_type b)
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

        vector_type negate(vector_type a) { //TODO }

        vector_type exp(vector_type a)
        {
            return _mm_set_ps(std::exp((float)a[0]),
                              std::exp((float)a[1]),
                              std::exp((float)a[2]),
                              std::exp((float)a[3]));
        }

        vector_type sqrt(vector_type a) { //TODO }
    };

}

#endif
