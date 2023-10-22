#ifndef MARRAY_BLAS_HPP
#define MARRAY_BLAS_HPP

#include "detail/fortran.h"

#ifdef MARRAY_USE_BLIS
#define BLIS_DISABLE_BLAS
#define _DEFINED_SCOMPLEX
#define _DEFINED_DCOMPLEX
#include <blis.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#define MARRAY_FOR_EACH_REAL_TYPE \
MARRAY_FOR_EACH_TYPE_BODY(float, float, float, S, s, S, s) \
MARRAY_FOR_EACH_TYPE_BODY(double, double, double, D, d, D, d)

#define MARRAY_FOR_EACH_COMPLEX_TYPE \
MARRAY_FOR_EACH_TYPE_BODY(scomplex, scomplex_f, float, C, c, S, s) \
MARRAY_FOR_EACH_TYPE_BODY(dcomplex, dcomplex_f, double, Z, z, D, d)

#define MARRAY_FOR_EACH_TYPE \
MARRAY_FOR_EACH_REAL_TYPE \
MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 1 BLAS, FORTRAN prototypes
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void    MARRAY_FC_FUNC(ch##rotg,CH##ROTG)      (ctype* a, ctype* b, ctyper* c, ctype* s); \
void    MARRAY_FC_FUNC(ch##swap,CH##SWAP)      (const integer* n,                           ctype* x, const integer* incx, ctype* y, const integer* incy); \
void    MARRAY_FC_FUNC(ch##scal,CH##SCAL)      (const integer* n, const ctype* alpha,       ctype* x, const integer* incx); \
void    MARRAY_FC_FUNC(ch##copy,CH##COPY)      (const integer* n,                     const ctype* x, const integer* incx, ctype* y, const integer* incy); \
void    MARRAY_FC_FUNC(ch##axpy,CH##AXPY)      (const integer* n, const ctype* alpha, const ctype* x, const integer* incx, ctype* y, const integer* incy); \
integer MARRAY_FC_FUNC(i##ch##amax,I##CH##AMAX)(const integer* n,                     const ctype* x, const integer* incx);

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void    MARRAY_FC_FUNC(ch##rotmg,CH##ROTMG)(ctype* d1, ctype* d2, ctype* a, const ctype* b, ctype* param); \
void    MARRAY_FC_FUNC(ch##rot,CH##ROT)    (const integer* n,       ctype* x, const integer* incx,       ctype* y, const integer* incy, const ctype* c, const ctype* s); \
void    MARRAY_FC_FUNC(ch##rotm,CH##ROTM)  (const integer* n,       ctype* x, const integer* incx,       ctype* y, const integer* incy, ctype* param); \
ctypef  MARRAY_FC_FUNC(ch##dot,CH##DOT)    (const integer* n, const ctype* x, const integer* incx, const ctype* y, const integer* incy); \
ctypef  MARRAY_FC_FUNC(ch##nrm2,CH##NRM2)  (const integer* n, const ctype* x, const integer* incx); \
ctypef  MARRAY_FC_FUNC(ch##asum,CH##ASUM)  (const integer* n, const ctype* x, const integer* incx);

MARRAY_FOR_EACH_REAL_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void   MARRAY_FC_FUNC(ch##chr##rot,CH##CHR##ROT)  (const integer* n,                            ctype* x, const integer* incx,       ctype* y, const integer* incy, const ctyper* c, const ctyper* s); \
void   MARRAY_FC_FUNC(ch##chr##scal,CH##CHR##SCAL)(const integer* n, const ctyper* alpha,       ctype* x, const integer* incx); \
ctypef MARRAY_FC_FUNC(ch##dotu,CH##DOTU)          (const integer* n,                      const ctype* x, const integer* incx, const ctype* y, const integer* incy); \
ctypef MARRAY_FC_FUNC(ch##dotc,CH##DOTC)          (const integer* n,                      const ctype* x, const integer* incx, const ctype* y, const integer* incy); \
ctyper MARRAY_FC_FUNC(chr##ch##nrm2,CHR##CH##NRM2)(const integer* n,                      const ctype* x, const integer* incx); \
ctyper MARRAY_FC_FUNC(chr##ch##asum,CHR##CH##ASUM)(const integer* n,                      const ctype* x, const integer* incx);

MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 2 BLAS, FORTRAN prototypes
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void MARRAY_FC_FUNC(ch##gemv,CH##GEMV)(                  const char* trans,                   const integer* m, const integer* n,                                       const ctype* alpha, const ctype* a, const integer* lda,  const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##gbmv,CH##GBMV)(                  const char* trans,                   const integer* m, const integer* n, const integer* kl, const integer* ku, const ctype* alpha, const ctype* a, const integer* lda,  const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##trmv,CH##TRMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const ctype* a, const integer* lda,        ctype* x, const integer* incx); \
void MARRAY_FC_FUNC(ch##tbmv,CH##TBMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                         const ctype* a, const integer* lda,        ctype* x, const integer* incx); \
void MARRAY_FC_FUNC(ch##tpmv,CH##TPMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const ctype* ap,                           ctype* x, const integer* incx); \
void MARRAY_FC_FUNC(ch##trsv,CH##TRSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const ctype* a, const integer* lda,        ctype* x, const integer* incx); \
void MARRAY_FC_FUNC(ch##tbsv,CH##TBSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                         const ctype* a, const integer* lda,        ctype* x, const integer* incx); \
void MARRAY_FC_FUNC(ch##tpsv,CH##TPSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const ctype* ap,                           ctype* x, const integer* incx); \

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void MARRAY_FC_FUNC(ch##symv,CH##SYMV)(const char* uplo,                                                        const integer* n,                                       const ctype* alpha, const ctype* a, const integer* lda,  const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##sbmv,CH##SBMV)(const char* uplo,                                                        const integer* n, const integer* k,                     const ctype* alpha, const ctype* a, const integer* lda,  const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##spmv,CH##SPMV)(const char* uplo,                                                        const integer* n,                                       const ctype* alpha, const ctype* ap,                     const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##ger,CH##GER)  (                                                       const integer* m, const integer* n,                                       const ctype* alpha, const ctype* x, const integer* incx, const ctype* y, const integer* incy, ctype* a, const integer* lda); \
void MARRAY_FC_FUNC(ch##syr,CH##SYR)  (const char* uplo,                                                        const integer* n,                                       const ctype* alpha, const ctype* x, const integer* incx,                                      ctype* a, const integer* lda); \
void MARRAY_FC_FUNC(ch##spr,CH##SPR)  (const char* uplo,                                                        const integer* n,                                       const ctype* alpha, const ctype* x, const integer* incx,                                      ctype* ap); \
void MARRAY_FC_FUNC(ch##syr2,CH##SYR2)(const char* uplo,                                                        const integer* n,                                       const ctype* alpha, const ctype* x, const integer* incx, const ctype* y, const integer* incy, ctype* a, const integer* lda); \
void MARRAY_FC_FUNC(ch##spr2,CH##SPR2)(const char* uplo,                                                        const integer* n,                                       const ctype* alpha, const ctype* x, const integer* incx, const ctype* y, const integer* incy, ctype* ap);

MARRAY_FOR_EACH_REAL_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void MARRAY_FC_FUNC(ch##hemv,CH##HEMV)(const char* uplo,                                                        const integer* n,                                       const ctype*  alpha, const ctype* a, const integer* lda,  const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##hbmv,CH##HBMV)(const char* uplo,                                                        const integer* n, const integer* k,                     const ctype*  alpha, const ctype* a, const integer* lda,  const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##hpmv,CH##HPMV)(const char* uplo,                                                        const integer* n,                                       const ctype*  alpha, const ctype* ap,                     const ctype* x, const integer* incx, const ctype* beta, ctype* y, const integer* incy); \
void MARRAY_FC_FUNC(ch##geru,CH##GERU)(                                                       const integer* m, const integer* n,                                       const ctype*  alpha, const ctype* x, const integer* incx, const ctype* y, const integer* incy, ctype* a, const integer* lda); \
void MARRAY_FC_FUNC(ch##gerc,CH##GERC)(                                                       const integer* m, const integer* n,                                       const ctype*  alpha, const ctype* x, const integer* incx, const ctype* y, const integer* incy, ctype* a, const integer* lda); \
void MARRAY_FC_FUNC(ch##her,CH##HER)  (const char* uplo,                                                        const integer* n,                                       const ctyper* alpha, const ctype* x, const integer* incx,                                      ctype* a, const integer* lda); \
void MARRAY_FC_FUNC(ch##hpr,CH##HPR)  (const char* uplo,                                                        const integer* n,                                       const ctyper* alpha, const ctype* x, const integer* incx,                                      ctype* ap); \
void MARRAY_FC_FUNC(ch##her2,CH##HER2)(const char* uplo,                                                        const integer* n,                                       const ctype*  alpha, const ctype* x, const integer* incx, const ctype* y, const integer* incy, ctype* a, const integer* lda); \
void MARRAY_FC_FUNC(ch##hpr2,CH##HPR2)(const char* uplo,                                                        const integer* n,                                       const ctype*  alpha, const ctype* x, const integer* incx, const ctype* y, const integer* incy, ctype* ap);

MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 3 BLAS, FORTRAN prototypes
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void MARRAY_FC_FUNC(ch##gemm,CH##GEMM)  (                                    const char* transa, const char* transb,                   const integer* m, const integer* n, const integer* k, const ctype* alpha, const ctype* a, const integer* lda, const ctype* b, const integer* ldb, const ctype* beta, ctype* c, const integer* ldc); \
void MARRAY_FC_FUNC(ch##symm,CH##SYMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const ctype* alpha, const ctype* a, const integer* lda, const ctype* b, const integer* ldb, const ctype* beta, ctype* c, const integer* ldc); \
void MARRAY_FC_FUNC(ch##syrk,CH##SYRK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const ctype* alpha, const ctype* a, const integer* lda,                                     const ctype* beta, ctype* c, const integer* ldc); \
void MARRAY_FC_FUNC(ch##syr2k,CH##SYR2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const ctype* alpha, const ctype* a, const integer* lda, const ctype* b, const integer* ldb, const ctype* beta, ctype* c, const integer* ldc); \
void MARRAY_FC_FUNC(ch##trmm,CH##TRMM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const ctype* alpha, const ctype* a, const integer* lda,       ctype* b, const integer* ldb); \
void MARRAY_FC_FUNC(ch##trsm,CH##TRSM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const ctype* alpha, const ctype* a, const integer* lda,       ctype* b, const integer* ldb);

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
void MARRAY_FC_FUNC(ch##hemm,CH##HEMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const ctype*  alpha, const ctype* a, const integer* lda, const ctype* b, const integer* ldb, const ctype*  beta, ctype* c, const integer* ldc); \
void MARRAY_FC_FUNC(ch##herk,CH##HERK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const ctyper* alpha, const ctype* a, const integer* lda,                                     const ctyper* beta, ctype* c, const integer* ldc); \
void MARRAY_FC_FUNC(ch##her2k,CH##HER2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const ctype*  alpha, const ctype* a, const integer* lda, const ctype* b, const integer* ldb, const ctyper* beta, ctype* c, const integer* ldc);

MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 1 BLAS, C wrappers
 *
 *****************************************************************************/

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline void   c_##ch##rotg (ctype* a, ctype* b, ctyper* c, ctype* s) \
{ \
    MARRAY_FC_FUNC(ch##rotg,CH##ROTG)(a, b, c, s); \
} \
\
static inline void   c_##ch##swap (const integer n, ctype* x, const integer incx, \
                                                    ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##swap,CH##SWAP)(&n, x, &incx, y, &incy); \
} \
\
static inline void   c_##ch##scal (const integer n, const ctype alpha, ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##scal,CH##SCAL)(&n, &alpha, x, &incx); \
} \
\
static inline void   c_##ch##copy (const integer n, const ctype* x, const integer incx, \
                                                          ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##copy,CH##COPY)(&n, x, &incx, y, &incy); \
} \
\
static inline void   c_##ch##axpy (const integer n, const ctype alpha, const ctype* x, const integer incx, \
                                                                             ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##axpy,CH##AXPY)(&n, &alpha, x, &incx, y, &incy); \
} \
\
static inline integer c_i##ch##amax(const integer n, const ctype* x, const integer incx) \
{ \
    return MARRAY_FC_FUNC(i##ch##amax,I##CH##AMAX)(&n, x, &incx)-1; \
}

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline ctype c_##ch##dot  (const integer n, const ctype* x, const integer incx, \
                                              const ctype* y, const integer incy) \
{ \
    return MARRAY_FC_FUNC(ch##dot,CH##DOT)(&n, x, &incx, y, &incy); \
} \
\
static inline ctype c_##ch##nrm2 (const integer n, const ctype* x, const integer incx) \
{ \
    return MARRAY_FC_FUNC(ch##nrm2,CH##NRM2)(&n, x, &incx); \
} \
\
static inline ctype c_##ch##asum (const integer n, const ctype* x, const integer incx) \
{ \
    return MARRAY_FC_FUNC(ch##asum,CH##ASUM)(&n, x, &incx); \
} \
\
static inline void   c_##ch##rotmg(ctype* d1, ctype* d2, ctype* a, const ctype b, ctype* param) \
{ \
    MARRAY_FC_FUNC(ch##rotmg,CH##ROTMG)(d1, d2, a, &b, param); \
} \
\
static inline void   c_##ch##rot  (const integer n, ctype* x, const integer incx, \
                                               ctype* y, const integer incy, const ctype c, const ctype s) \
{ \
    MARRAY_FC_FUNC(ch##rot,CH##ROT)(&n, x, &incx, y, &incy, &c, &s); \
} \
\
static inline void   c_##ch##rotm (const integer n, ctype* x, const integer incx, \
                                               ctype* y, const integer incy, ctype* param) \
{ \
    MARRAY_FC_FUNC(ch##rotm,CH##ROTM)(&n, x, &incx, y, &incy, param); \
}

MARRAY_FOR_EACH_REAL_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline void   c_##ch##chr##rot  (const integer n, ctype* x, const integer incx, \
                                                         ctype* y, const integer incy, const ctyper c, const ctyper s) \
{ \
    MARRAY_FC_FUNC(ch##chr##rot,CH##CHR##ROT)(&n, x, &incx, y, &incy, &c, &s); \
} \
\
static inline void   c_##ch##chr##scal (const integer n, const ctyper alpha, ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##chr##scal,CH##CHR##SCAL)(&n, &alpha, x, &incx); \
} \
\
static inline ctype c_##ch##dotu (const integer n, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy) \
{ \
    ctypef tmp = MARRAY_FC_FUNC(ch##dotu,CH##DOTU)(&n, x, &incx, y, &incy); \
    return MAKE_COMPLEX(ctyper, tmp.real, tmp.imag); \
} \
\
static inline ctype c_##ch##dotc (const integer n, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy) \
{ \
    ctypef tmp = MARRAY_FC_FUNC(ch##dotc,CH##DOTC)(&n, x, &incx, y, &incy); \
    return MAKE_COMPLEX(ctyper, tmp.real, tmp.imag); \
} \
\
static inline ctyper c_##chr##ch##nrm2 (const integer n, const ctype* x, const integer incx) \
{ \
    return MARRAY_FC_FUNC(chr##ch##nrm2,CHR##CH##NRM2)(&n, x, &incx); \
} \
\
static inline ctyper c_##chr##ch##asum (const integer n, const ctype* x, const integer incx) \
{ \
    return MARRAY_FC_FUNC(chr##ch##asum,CHR##CH##ASUM)(&n, x, &incx); \
}

MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 2 BLAS, C wrappers
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline void c_##ch##gemv(const char trans, const integer m, const integer n, \
                                const ctype alpha, const ctype* a, const integer lda, \
                                                   const ctype* x, const integer incx, \
                                const ctype  beta,       ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##gemv,ZGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##gbmv(const char trans, \
                                const integer m, const integer n, const integer kl, const integer ku, \
                                const ctype alpha, const ctype* a, const integer lda, \
                                                   const ctype* x, const integer incx, \
                                const ctype  beta,       ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##gbmv,ZGBMV)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##trmv(const char uplo, const char trans, const char diag, const integer n, \
                                const ctype* a, const integer lda, \
                                      ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##trmv,ZTRMV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); \
} \
\
static inline void c_##ch##tbmv(const char uplo, const char trans, const char diag, \
                                const integer n, const integer k, \
                                const ctype* a, const integer lda, \
                                      ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##tbmv,ZTBMV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); \
} \
\
static inline void c_##ch##tpmv(const char uplo, const char trans, const char diag, const integer n, \
                                const ctype* ap, ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##tpmv,ZTPMV)(&uplo, &trans, &diag, &n, ap, x, &incx); \
} \
\
static inline void c_##ch##trsv(const char uplo, const char trans, const char diag, const integer n, \
                                const ctype* a, const integer lda, \
                                      ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##trsv,ZTRSV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); \
} \
\
static inline void c_##ch##tbsv(const char uplo, const char trans, const char diag, \
                                const integer n, const integer k, \
                                const ctype* a, const integer lda, \
                                      ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##tbsv,ZTBSV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); \
} \
\
static inline void c_##ch##tpsv(const char uplo, const char trans, const char diag, const integer n, \
                                const ctype* ap, ctype* x, const integer incx) \
{ \
    MARRAY_FC_FUNC(ch##tpsv,ZTPSV)(&uplo, &trans, &diag, &n, ap, x, &incx); \
}

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline void c_##ch##symv(const char uplo, const integer n, \
                                const ctype alpha, const ctype* a, const integer lda, \
                                                   const ctype* x, const integer incx, \
                                const ctype  beta,       ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##symv,SSYMV)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##sbmv(const char uplo, const integer n, const integer k, \
                                const ctype alpha, const ctype* a, const integer lda, \
                                                   const ctype* x, const integer incx, \
                                const ctype  beta,       ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##sbmv,SSBMV)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##spmv(const char uplo, const integer n, \
                                const ctype alpha, const ctype* ap, \
                                                   const ctype*  x, const integer incx, \
                                const ctype  beta,       ctype*  y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##spmv,SSPMV)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##ger (const integer m, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy, \
                                                         ctype* a, const integer lda) \
{ \
    MARRAY_FC_FUNC(ch##ger,SGER)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); \
} \
\
static inline void c_##ch##syr (const char uplo, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                         ctype* a, const integer lda) \
{ \
    MARRAY_FC_FUNC(ch##syr,CH##SYR)(&uplo, &n, &alpha, x, &incx, a, &lda); \
} \
\
static inline void c_##ch##spr (const char uplo, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                         ctype* ap) \
{ \
    MARRAY_FC_FUNC(ch##spr,CH##SPR)(&uplo, &n, &alpha, x, &incx, ap); \
} \
\
static inline void c_##ch##syr2(const char uplo, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy, \
                                                         ctype* a, const integer lda) \
{ \
    MARRAY_FC_FUNC(ch##syr2,CH##SYR2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda); \
} \
\
static inline void c_##ch##spr2(const char uplo, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy, \
                                                         ctype* ap) \
{ \
    MARRAY_FC_FUNC(ch##spr2,CH##SPR2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap); \
}

MARRAY_FOR_EACH_REAL_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline void c_##ch##hemv(const char uplo, const integer n, \
                                const ctype alpha, const ctype* a, const integer lda, \
                                                   const ctype* x, const integer incx, \
                                const ctype  beta,       ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##hemv,CH##HEMV)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##hbmv(const char uplo, const integer n, const integer k, \
                                const ctype alpha, const ctype* a, const integer lda, \
                                                   const ctype* x, const integer incx, \
                                const ctype  beta,       ctype* y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##hbmv,CH##HBMV)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##hpmv(const char uplo, const integer n, \
                                const ctype alpha, const ctype* ap, \
                                                   const ctype*  x, const integer incx, \
                                const ctype  beta,       ctype*  y, const integer incy) \
{ \
    MARRAY_FC_FUNC(ch##hpmv,CH##HPMV)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); \
} \
\
static inline void c_##ch##geru(const integer m, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy, \
                                                         ctype* a, const integer lda) \
{ \
    MARRAY_FC_FUNC(ch##geru,CH##GERU)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); \
} \
\
static inline void c_##ch##gerc(const integer m, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy, \
                                                         ctype* a, const integer lda) \
{ \
    MARRAY_FC_FUNC(ch##gerc,CH##GERC)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); \
} \
\
static inline void c_##ch##her (const char uplo, const integer n, \
                                const ctyper alpha, const ctype* x, const integer incx, \
                                                          ctype* a, const integer lda) \
{ \
    MARRAY_FC_FUNC(ch##her,CH##HER)(&uplo, &n, &alpha, x, &incx, a, &lda); \
} \
\
static inline void c_##ch##hpr (const char uplo, const integer n, \
                                const ctyper alpha, const ctype* x, const integer incx, \
                                                          ctype* ap) \
{ \
    MARRAY_FC_FUNC(ch##hpr,CH##HPR)(&uplo, &n, &alpha, x, &incx, ap); \
} \
\
static inline void c_##ch##her2(const char uplo, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy, \
                                                         ctype* a, const integer lda) \
{ \
    MARRAY_FC_FUNC(ch##her2,CH##HER2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda); \
} \
\
static inline void c_##ch##hpr2(const char uplo, const integer n, \
                                const ctype alpha, const ctype* x, const integer incx, \
                                                   const ctype* y, const integer incy, \
                                                         ctype* ap) \
{ \
    MARRAY_FC_FUNC(ch##hpr2,CH##HPR2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap); \
}

MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 3 BLAS, C wrappers
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline void c_##ch##gemm (const char transa, const char transb, \
                                 const integer m, const integer n, const integer k, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                                    const ctype* b, const integer ldb, \
                                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    MARRAY_FC_FUNC(ch##gemm,CH##GEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
} \
\
static inline void c_##ch##symm (const char side, const char uplo, \
                                 const integer m, const integer n, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                                    const ctype* b, const integer ldb, \
                                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    MARRAY_FC_FUNC(ch##symm,CH##SYMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
} \
\
static inline void c_##ch##syrk (const char uplo, const char trans, \
                                 const integer n, const integer k, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    MARRAY_FC_FUNC(ch##syrk,CH##SYRK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); \
} \
\
static inline void c_##ch##syr2k(const char uplo, const char trans, \
                                 const integer n, const integer k, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                                    const ctype* b, const integer ldb, \
                                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    MARRAY_FC_FUNC(ch##syr2k,CH##SYR2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
} \
\
static inline void c_##ch##trmm (const char side, const char uplo, const char transa, const char diag, \
                                 const integer m, const integer n, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                                          ctype* b, const integer ldb) \
{ \
    MARRAY_FC_FUNC(ch##trmm,CH##TRMM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); \
} \
\
static inline void c_##ch##trsm (const char side, const char uplo, const char transa, const char diag, \
                                 const integer m, const integer n, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                                          ctype* b, const integer ldb) \
{ \
    MARRAY_FC_FUNC(ch##trsm,CH##TRSM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); \
}

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
static inline void c_##ch##hemm (const char side, const char uplo, \
                                 const integer m, const integer n, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                                    const ctype* b, const integer ldb, \
                                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    MARRAY_FC_FUNC(ch##hemm,CH##HEMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
} \
\
static inline void c_##ch##herk (const char uplo, const char trans, \
                                 const integer n, const integer k, \
                                 const ctyper alpha, const ctype* a, const integer lda, \
                                 const ctyper  beta,       ctype* c, const integer ldc) \
{ \
    MARRAY_FC_FUNC(ch##herk,CH##HERK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); \
} \
\
static inline void c_##ch##her2k(const char uplo, const char trans, \
                                 const integer n, const integer k, \
                                 const ctype alpha, const ctype* a, const integer lda, \
                                                    const ctype* b, const integer ldb, \
                                 const ctyper beta,       ctype* c, const integer ldc) \
{ \
    MARRAY_FC_FUNC(ch##her2k,CH##HER2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
}

MARRAY_FOR_EACH_COMPLEX_TYPE

#ifdef __cplusplus
}
#endif

/*
 * #define more familiar names for the C versions
 */
#define srotg  c_srotg
#define srotmg c_srotmg
#define srot   c_srot
#define srotm  c_srotm
#define sswap  c_sswap
#define sscal  c_sscal
#define scopy  c_scopy
#define saxpy  c_saxpy
#define sdot   c_sdot
#define snrm2  c_snrm2
#define sasum  c_sasum
#define isamax c_isamax
#define sgemv  c_sgemv
#define sgbmv  c_sgbmv
#define ssymv  c_ssymv
#define ssbmv  c_ssbmv
#define sspmv  c_sspmv
#define strmv  c_strmv
#define stbmv  c_stbmv
#define stpmv  c_stpmv
#define strsv  c_strsv
#define stbsv  c_stbsv
#define stpsv  c_stpsv
#define sger   c_sger
#define ssyr   c_ssyr
#define sspr   c_sspr
#define ssyr2  c_ssyr2
#define sspr2  c_sspr2
#define sgemm  c_sgemm
#define ssymm  c_ssymm
#define ssyrk  c_ssyrk
#define ssyr2k c_ssyr2k
#define strmm  c_strmm
#define strsm  c_strsm

#define drotg  c_drotg
#define drotmg c_drotmg
#define drot   c_drot
#define drotm  c_drotm
#define dswap  c_dswap
#define dscal  c_dscal
#define dcopy  c_dcopy
#define daxpy  c_daxpy
#define ddot   c_ddot
#define dnrm2  c_dnrm2
#define dasum  c_dasum
#define idamax c_idamax
#define dgemv  c_dgemv
#define dgbmv  c_dgbmv
#define dsymv  c_dsymv
#define dsbmv  c_dsbmv
#define dspmv  c_dspmv
#define dtrmv  c_dtrmv
#define dtbmv  c_dtbmv
#define dtpmv  c_dtpmv
#define dtrsv  c_dtrsv
#define dtbsv  c_dtbsv
#define dtpsv  c_dtpsv
#define dger   c_dger
#define dsyr   c_dsyr
#define dspr   c_dspr
#define dsyr2  c_dsyr2
#define dspr2  c_dspr2
#define dgemm  c_dgemm
#define dsymm  c_dsymm
#define dsyrk  c_dsyrk
#define dsyr2k c_dsyr2k
#define dtrmm  c_dtrmm
#define dtrsm  c_dtrsm

#define crotg  c_crotg
#define csrot  c_csrot
#define cswap  c_cswap
#define cscal  c_cscal
#define csscal c_csscal
#define ccopy  c_ccopy
#define caxpy  c_caxpy
#define cdotu  c_cdotu
#define cdotc  c_cdotc
#define scnrm2 c_scnrm2
#define scasum c_scasum
#define icamax c_icamax
#define cgemv  c_cgemv
#define cgbmv  c_cgbmv
#define chemv  c_chemv
#define chbmv  c_chbmv
#define chpmv  c_chpmv
#define ctrmv  c_ctrmv
#define ctbmv  c_ctbmv
#define ctpmv  c_ctpmv
#define ctrsv  c_ctrsv
#define ctbsv  c_ctbsv
#define ctpsv  c_ctpsv
#define cgeru  c_cgeru
#define cgerc  c_cgerc
#define cher   c_cher
#define chpr   c_chpr
#define cher2  c_cher2
#define chpr2  c_chpr2
#define cgemm  c_cgemm
#define csymm  c_csymm
#define chemm  c_chemm
#define csyrk  c_csyrk
#define csyr2k c_csyr2k
#define cherk  c_cherk
#define cher2k c_cher2k
#define ctrmm  c_ctrmm
#define ctrsm  c_ctrsm

#define zrotg  c_zrotg
#define zdrot  c_zdrot
#define zswap  c_zswap
#define zscal  c_zscal
#define zdscal c_zdscal
#define zcopy  c_zcopy
#define zaxpy  c_zaxpy
#define zdotu  c_zdotu
#define zdotc  c_zdotc
#define dznrm2 c_dznrm2
#define dzasum c_dzasum
#define izamax c_izamax
#define zgemv  c_zgemv
#define zgbmv  c_zgbmv
#define zhemv  c_zhemv
#define zhbmv  c_zhbmv
#define zhpmv  c_zhpmv
#define ztrmv  c_ztrmv
#define ztbmv  c_ztbmv
#define ztpmv  c_ztpmv
#define ztrsv  c_ztrsv
#define ztbsv  c_ztbsv
#define ztpsv  c_ztpsv
#define zgeru  c_zgeru
#define zgerc  c_zgerc
#define zher   c_zher
#define zhpr   c_zhpr
#define zher2  c_zher2
#define zhpr2  c_zhpr2
#define zgemm  c_zgemm
#define zsymm  c_zsymm
#define zhemm  c_zhemm
#define zsyrk  c_zsyrk
#define zsyr2k c_zsyr2k
#define zherk  c_zherk
#define zher2k c_zher2k
#define ztrmm  c_ztrmm
#define ztrsm  c_ztrsm

#ifdef __cplusplus

#include "marray.hpp"

namespace MArray
{

namespace detail
{

template <typename T>
using value_type = std::remove_cv_t<typename std::decay_t<T>::value_type>;

}

namespace blas
{

/******************************************************************************
 *
 * Level 1 BLAS, C++ overloads
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void rotg(ctype* a, ctype* b, ctyper* c, ctype* s) \
{ \
    ch##rotg(a, b, c, s); \
} \
\
inline void swap(const integer n, ctype* x, const integer incx, \
                                  ctype* y, const integer incy) \
{ \
    ch##swap(n, x, incx, y, incy); \
} \
\
inline void scal(const integer n, const ctype alpha, ctype* x, const integer incx) \
{ \
    ch##scal(n, alpha, x, incx); \
} \
\
inline void copy(const integer n, const ctype* x, const integer incx, \
                                        ctype* y, const integer incy) \
{ \
    ch##copy(n, x, incx, y, incy); \
} \
\
inline void axpy(const integer n, const ctype alpha, const ctype* x, const integer incx, \
                                                           ctype* y, const integer incy) \
{ \
    ch##axpy(n, alpha, x, incx, y, incy); \
} \
\
inline integer iamax(const integer n, const ctype* x, const integer incx) \
{ \
    return i##ch##amax(n, x, incx); \
} \
\
inline ctyper amax(const integer n, const ctype* x, const integer incx) \
{ \
    return std::abs(x[i##ch##amax(n, x, incx)]); \
}

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void rotmg(ctype* d1, ctype* d2, ctype* a, const ctype b, ctype* param) \
{ \
    ch##rotmg(d1, d2, a, b, param); \
} \
\
inline void rot(const integer n, ctype* x, const integer incx, \
                                 ctype* y, const integer incy, const ctype c, const ctype s) \
{ \
    ch##rot(n, x, incx, y, incy, c, s); \
} \
\
inline void rotm(const integer n, ctype* x, const integer incx, \
                                  ctype* y, const integer incy, ctype* param) \
{ \
    ch##rotm(n, x, incx, y, incy, param); \
} \
\
inline ctype dot(const integer n, const ctype* x, const integer incx, \
                                  const ctype* y, const integer incy) \
{ \
    return ch##dot(n, x, incx, y, incy); \
} \
\
inline ctype dotc(const integer n, const ctype* x, const integer incx, \
                                   const ctype* y, const integer incy) \
{ \
    return ch##dot(n, x, incx, y, incy); \
} \
\
inline ctype dotu(const integer n, const ctype* x, const integer incx, \
                                   const ctype* y, const integer incy) \
{ \
    return ch##dot(n, x, incx, y, incy); \
} \
\
inline ctype nrm2(const integer n, const ctype* x, const integer incx) \
{ \
    return ch##nrm2(n, x, incx); \
} \
\
inline ctype asum(const integer n, const ctype* x, const integer incx) \
{ \
    return ch##asum(n, x, incx); \
}

MARRAY_FOR_EACH_REAL_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void rot(const integer n, ctype* x, const integer incx, \
                                 ctype* y, const integer incy, const ctyper c, const ctyper s) \
{ \
    ch##chr##rot(n, x, incx, y, incy, c, s); \
} \
\
inline void scal(const integer n, const ctyper alpha, ctype* x, const integer incx) \
{ \
    ch##chr##scal(n, alpha, x, incx); \
} \
\
inline ctype dot(const integer n, const ctype* x, const integer incx, \
                                  const ctype* y, const integer incy) \
{ \
    return ch##dotc(n, x, incx, y, incy); \
} \
\
inline ctype dotc(const integer n, const ctype* x, const integer incx, \
                                   const ctype* y, const integer incy) \
{ \
    return ch##dotc(n, x, incx, y, incy); \
} \
\
inline ctype dotu(const integer n, const ctype* x, const integer incx, \
                                   const ctype* y, const integer incy) \
{ \
    return ch##dotu(n, x, incx, y, incy); \
} \
\
inline ctyper nrm2(const integer n, const ctype* x, const integer incx) \
{ \
    return chr##ch##nrm2(n, x, incx); \
} \
\
inline ctyper asum(const integer n, const ctype* x, const integer incx) \
{ \
    return chr##ch##asum(n, x, incx); \
}

MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 2 BLAS, C++ overloads
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void gemv(const char trans, const integer m, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* x, const integer incx, \
                 const ctype  beta,       ctype* y, const integer incy) \
{ \
    ch##gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); \
} \
\
inline void gbmv(const char trans, \
                 const integer m, const integer n, const integer kl, const integer ku, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* x, const integer incx, \
                 const ctype  beta,       ctype* y, const integer incy) \
{ \
    ch##gbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy); \
} \
\
inline void trmv(const char uplo, const char trans, const char diag, const integer n, \
                 const ctype* a, const integer lda, \
                       ctype* x, const integer incx) \
{ \
    ch##trmv(uplo, trans, diag, n, a, lda, x, incx); \
} \
\
inline void tbmv(const char uplo, const char trans, const char diag, \
                 const integer n, const integer k, \
                 const ctype* a, const integer lda, \
                       ctype* x, const integer incx) \
{ \
    ch##tbmv(uplo, trans, diag, n, k, a, lda, x, incx); \
} \
\
inline void tpmv(const char uplo, const char trans, const char diag, const integer n, \
                 const ctype* ap, ctype* x, const integer incx) \
{ \
    ch##tpmv(uplo, trans, diag, n, ap, x, incx); \
} \
\
inline void trsv(const char uplo, const char trans, const char diag, const integer n, \
                 const ctype* a, const integer lda, \
                       ctype* x, const integer incx) \
{ \
    ch##trsv(uplo, trans, diag, n, a, lda, x, incx); \
} \
\
inline void tbsv(const char uplo, const char trans, const char diag, \
                 const integer n, const integer k, \
                 const ctype* a, const integer lda, \
                       ctype* x, const integer incx) \
{ \
    ch##tbsv(uplo, trans, diag, n, k, a, lda, x, incx); \
} \
\
inline void tpsv(const char uplo, const char trans, const char diag, const integer n, \
                 const ctype* ap, ctype* x, const integer incx) \
{ \
    ch##tpsv(uplo, trans, diag, n, ap, x, incx); \
}

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void hemv(const char uplo, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* x, const integer incx, \
                 const ctype  beta,       ctype* y, const integer incy) \
{ \
    ch##symv(uplo, n, alpha, a, lda, x, incx, beta, y, incy); \
} \
\
inline void hbmv(const char uplo, const integer n, const integer k, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* x, const integer incx, \
                 const ctype  beta,       ctype* y, const integer incy) \
{ \
    ch##sbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy); \
} \
\
inline void hpmv(const char uplo, const integer n, \
                 const ctype alpha, const ctype* ap, \
                                    const ctype*  x, const integer incx, \
                 const ctype  beta,       ctype*  y, const integer incy) \
{ \
    ch##spmv(uplo, n, alpha, ap, x, incx, beta, y, incy); \
} \
\
inline void ger(const integer m, const integer n, \
                const ctype alpha, const ctype* x, const integer incx, \
                                   const ctype* y, const integer incy, \
                                         ctype* a, const integer lda) \
{ \
    ch##ger(m, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void gerc(const integer m, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* a, const integer lda) \
{ \
    ch##ger(m, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void geru(const integer m, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* a, const integer lda) \
{ \
    ch##ger(m, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void her(const char uplo, const integer n, \
                const ctype alpha, const ctype* x, const integer incx, \
                                         ctype* a, const integer lda) \
{ \
    ch##syr(uplo, n, alpha, x, incx, a, lda); \
} \
\
inline void hpr(const char uplo, const integer n, \
                const ctype alpha, const ctype* x, const integer incx, \
                                         ctype* ap) \
{ \
    ch##spr(uplo, n, alpha, x, incx, ap); \
} \
\
inline void her2(const char uplo, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* a, const integer lda) \
{ \
    ch##syr2(uplo, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void hpr2(const char uplo, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* ap) \
{ \
    ch##spr2(uplo, n, alpha, x, incx, y, incy, ap); \
}

MARRAY_FOR_EACH_REAL_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void hemv(const char uplo, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* x, const integer incx, \
                 const ctype  beta,       ctype* y, const integer incy) \
{ \
    ch##hemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy); \
} \
\
inline void hbmv(const char uplo, const integer n, const integer k, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* x, const integer incx, \
                 const ctype  beta,       ctype* y, const integer incy) \
{ \
    ch##hbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy); \
} \
\
inline void hpmv(const char uplo, const integer n, \
                 const ctype alpha, const ctype* ap, \
                                    const ctype*  x, const integer incx, \
                 const ctype  beta,       ctype*  y, const integer incy) \
{ \
    ch##hpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); \
} \
\
inline void ger(const integer m, const integer n, \
                const ctype alpha, const ctype* x, const integer incx, \
                                   const ctype* y, const integer incy, \
                                         ctype* a, const integer lda) \
{ \
    ch##gerc(m, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void gerc(const integer m, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* a, const integer lda) \
{ \
    ch##gerc(m, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void geru(const integer m, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* a, const integer lda) \
{ \
    ch##geru(m, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void her(const char uplo, const integer n, \
                const ctyper alpha, const ctype* x, const integer incx, \
                                          ctype* a, const integer lda) \
{ \
    ch##her(uplo, n, alpha, x, incx, a, lda); \
} \
\
inline void hpr(const char uplo, const integer n, \
                const ctyper alpha, const ctype* x, const integer incx, \
                                          ctype* ap) \
{ \
    ch##hpr(uplo, n, alpha, x, incx, ap); \
} \
\
inline void her2(const char uplo, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* a, const integer lda) \
{ \
    ch##her2(uplo, n, alpha, x, incx, y, incy, a, lda); \
} \
\
inline void hpr2(const char uplo, const integer n, \
                 const ctype alpha, const ctype* x, const integer incx, \
                                    const ctype* y, const integer incy, \
                                          ctype* ap) \
{ \
    ch##hpr2(uplo, n, alpha, x, incx, y, incy, ap); \
}

MARRAY_FOR_EACH_COMPLEX_TYPE

/******************************************************************************
 *
 * Level 3 BLAS, C++ overloads
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void gemm(const char transa, const char transb, \
                 const integer m, const integer n, const integer k, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* b, const integer ldb, \
                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
} \
\
inline void symm(const char side, const char uplo, \
                 const integer m, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* b, const integer ldb, \
                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); \
} \
\
inline void syrk(const char uplo, const char trans, \
                 const integer n, const integer k, \
                 const ctype alpha, const ctype* a, const integer lda, \
                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc); \
} \
\
inline void syr2k(const char uplo, const char trans, \
                  const integer n, const integer k, \
                  const ctype alpha, const ctype* a, const integer lda, \
                                     const ctype* b, const integer ldb, \
                  const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##syr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
} \
\
inline void trmm(const char side, const char uplo, const char transa, const char diag, \
                 const integer m, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                          ctype* b, const integer ldb) \
{ \
    ch##trmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb); \
} \
\
inline void trsm(const char side, const char uplo, const char transa, const char diag, \
                 const integer m, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                          ctype* b, const integer ldb) \
{ \
    ch##trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb); \
}

MARRAY_FOR_EACH_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void hemm(const char side, const char uplo, \
                 const integer m, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* b, const integer ldb, \
                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); \
} \
\
inline void herk(const char uplo, const char trans, \
                 const integer n, const integer k, \
                 const ctype alpha, const ctype* a, const integer lda, \
                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc); \
} \
\
inline void her2k(const char uplo, const char trans, \
                  const integer n, const integer k, \
                  const ctype alpha, const ctype* a, const integer lda, \
                                     const ctype* b, const integer ldb, \
                  const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##syr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
}

MARRAY_FOR_EACH_REAL_TYPE

#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void hemm(const char side, const char uplo, \
                 const integer m, const integer n, \
                 const ctype alpha, const ctype* a, const integer lda, \
                                    const ctype* b, const integer ldb, \
                 const ctype  beta,       ctype* c, const integer ldc) \
{ \
    ch##hemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); \
} \
\
inline void herk(const char uplo, const char trans, \
                 const integer n, const integer k, \
                 const ctyper alpha, const ctype* a, const integer lda, \
                 const ctyper  beta,       ctype* c, const integer ldc) \
{ \
    ch##herk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc); \
} \
\
inline void her2k(const char uplo, const char trans, \
                  const integer n, const integer k, \
                  const ctype alpha, const ctype* a, const integer lda, \
                                     const ctype* b, const integer ldb, \
                  const ctyper beta,       ctype* c, const integer ldc) \
{ \
    ch##her2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
}

MARRAY_FOR_EACH_COMPLEX_TYPE

#ifdef BLIS_H

/******************************************************************************
 *
 * Level 1 BLIS, C++ overloads
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void bli_addv \
      ( \
              conj_t conjx, \
              dim_t  n, \
        const ctype* x, inc_t incx, \
              ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##addv( conjx, n, x, incx, y, incy ); \
} \
\
inline void bli_copyv \
      ( \
              conj_t conjx, \
              dim_t  n, \
        const ctype* x, inc_t incx, \
              ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##copyv( conjx, n, x, incx, y, incy ); \
} \
\
inline void bli_subv \
      ( \
              conj_t conjx, \
              dim_t  n, \
        const ctype* x, inc_t incx, \
              ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##subv( conjx, n, x, incx, y, incy ); \
} \
\
inline dim_t bli_amaxv \
     ( \
             dim_t  n, \
       const ctype* x, inc_t incx \
      ) \
{ \
    dim_t index; \
    bli_##ch##amaxv( n, x, incx, &index ); \
    return index; \
} \
\
inline void bli_axpbyv \
     ( \
             conj_t conjx, \
             dim_t  n, \
       const ctype& alpha, \
       const ctype* x, inc_t incx, \
       const ctype& beta, \
             ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##axpbyv( conjx, n, &alpha, x, incx, &beta, y, incy ); \
} \
\
inline void bli_axpyv \
     ( \
             conj_t conjx, \
             dim_t  n, \
       const ctype& alpha, \
       const ctype* x, inc_t incx, \
             ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##axpyv( conjx, n, &alpha, x, incx, y, incy ); \
} \
\
inline void bli_scal2v \
     ( \
             conj_t conjx, \
             dim_t  n, \
       const ctype& alpha, \
       const ctype* x, inc_t incx, \
             ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##scal2v( conjx, n, &alpha, x, incx, y, incy ); \
} \
\
inline ctype bli_dotv \
     ( \
             conj_t conjx, \
             conj_t conjy, \
             dim_t  n, \
       const ctype* x, inc_t incx, \
       const ctype* y, inc_t incy \
      ) \
{ \
    ctype rho; \
    bli_##ch##dotv( conjx, conjy, n, x, incx, y, incy, &rho ); \
    return rho; \
} \
\
inline void bli_dotxv \
     ( \
             conj_t conjx, \
             conj_t conjy, \
             dim_t  n, \
       const ctype& alpha, \
       const ctype* x, inc_t incx, \
       const ctype* y, inc_t incy, \
       const ctype& beta, \
             ctype& rho  \
      ) \
{ \
    bli_##ch##dotxv( conjx, conjy, n, &alpha, x, incx, y, incy, &beta, &rho ); \
} \
\
inline void bli_invertv \
     ( \
       dim_t  n, \
       ctype* x, inc_t incx  \
      ) \
{ \
    bli_##ch##invertv( n, x, incx ); \
} \
\
inline void bli_invscalv \
     ( \
             conj_t conjalpha, \
             dim_t  n, \
       const ctype& alpha, \
             ctype* x, inc_t incx  \
      ) \
{ \
    bli_##ch##invscalv( conjalpha, n, &alpha, x, incx ); \
} \
\
inline void bli_scalv \
     ( \
             conj_t conjalpha, \
             dim_t  n, \
       const ctype& alpha, \
             ctype* x, inc_t incx  \
      ) \
{ \
    bli_##ch##scalv( conjalpha, n, &alpha, x, incx ); \
} \
\
inline void bli_setv \
     ( \
             conj_t conjalpha, \
             dim_t  n, \
       const ctype& alpha, \
             ctype* x, inc_t incx  \
      ) \
{ \
    bli_##ch##setv( conjalpha, n, &alpha, x, incx ); \
} \
\
inline void bli_swapv \
     ( \
       dim_t  n, \
       ctype* x, inc_t incx, \
       ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##swapv( n, x, incx, y, incy ); \
} \
\
inline void bli_xpbyv \
     ( \
             conj_t conjx, \
             dim_t  n, \
       const ctype* x, inc_t incx, \
       const ctype& beta, \
             ctype* y, inc_t incy  \
      ) \
{ \
    bli_##ch##xpbyv( conjx, n, x, incx, &beta, y, incy ); \
} \
\
inline ctyper bli_asumv \
     ( \
             dim_t    n, \
       const ctype*   x, inc_t incx \
     ) \
{ \
    ctyper asum; \
    bli_##ch##asumv( n, x, incx, &asum ); \
    return asum; \
} \
\
inline ctyper norm1v \
     ( \
             dim_t    n, \
       const ctype*   x, inc_t incx \
     ) \
{ \
    ctyper norm; \
    bli_##ch##norm1v( n, x, incx, &norm ); \
    return norm; \
} \
\
inline ctyper normfv \
     ( \
             dim_t    n, \
       const ctype*   x, inc_t incx \
     ) \
{ \
    ctyper norm; \
    bli_##ch##normfv( n, x, incx, &norm ); \
    return norm; \
} \
\
inline ctyper normiv \
     ( \
             dim_t    n, \
       const ctype*   x, inc_t incx \
     ) \
{ \
    ctyper norm; \
    bli_##ch##normiv( n, x, incx, &norm ); \
    return norm; \
}

MARRAY_FOR_EACH_TYPE

/******************************************************************************
 *
 * Level 2 BLAS, C++ overloads
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void bli_gemv \
     ( \
             trans_t transa, \
             conj_t  conjx, \
             dim_t   m, \
             dim_t   n, \
       const ctype*  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  x, inc_t incx, \
       const ctype*  beta, \
             ctype*  y, inc_t incy  \
     ) \
{ \
    bli_##ch##gemv( transa, conjx, m, n, alpha, a, rs_a, cs_a, x, incx, beta, y, incy ); \
} \
\
inline void bli_ger \
     ( \
             conj_t conjx, \
             conj_t conjy, \
             dim_t  m, \
             dim_t  n, \
       const ctype* alpha, \
       const ctype* x, inc_t incx, \
       const ctype* y, inc_t incy, \
             ctype* a, inc_t rs_a, inc_t cs_a  \
     ) \
{ \
    bli_##ch##ger( conjx, conjy, m, n, alpha, x, incx, y, incy, a, rs_a, cs_a ); \
} \
\
inline void bli_hemv \
     ( \
             uplo_t uploa, \
             conj_t conja, \
             conj_t conjx, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* a, inc_t rs_a, inc_t cs_a, \
       const ctype* x, inc_t incx, \
       const ctype* beta, \
             ctype* y, inc_t incy  \
     ) \
{ \
    bli_##ch##hemv( uploa, conja, conjx, m, alpha, a, rs_a, cs_a, x, incx, beta, y, incy ); \
} \
\
inline void bli_symv \
     ( \
             uplo_t uploa, \
             conj_t conja, \
             conj_t conjx, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* a, inc_t rs_a, inc_t cs_a, \
       const ctype* x, inc_t incx, \
       const ctype* beta, \
             ctype* y, inc_t incy  \
     ) \
{ \
    bli_##ch##symv( uploa, conja, conjx, m, alpha, a, rs_a, cs_a, x, incx, beta, y, incy ); \
} \
\
inline void bli_shmv \
     ( \
             uplo_t uploa, \
             conj_t conja, \
             conj_t conjx, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* a, inc_t rs_a, inc_t cs_a, \
       const ctype* x, inc_t incx, \
       const ctype* beta, \
             ctype* y, inc_t incy  \
     ) \
{ \
    bli_##ch##shmv( uploa, conja, conjx, m, alpha, a, rs_a, cs_a, x, incx, beta, y, incy ); \
} \
\
inline void bli_skmv \
     ( \
             uplo_t uploa, \
             conj_t conja, \
             conj_t conjx, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* a, inc_t rs_a, inc_t cs_a, \
       const ctype* x, inc_t incx, \
       const ctype* beta, \
             ctype* y, inc_t incy  \
     ) \
{ \
    bli_##ch##skmv( uploa, conja, conjx, m, alpha, a, rs_a, cs_a, x, incx, beta, y, incy ); \
} \
\
inline void bli_her \
     ( \
             uplo_t   uploa, \
             conj_t   conjx, \
             dim_t    m, \
       const ctyper* alpha, \
       const ctype*   x, inc_t incx, \
             ctype*   a, inc_t rs_a, inc_t cs_a  \
     ) \
{ \
    bli_##ch##her( uploa, conjx, m, alpha, x, incx, a, rs_a, cs_a ); \
} \
\
inline void bli_syr \
     ( \
             uplo_t uploa, \
             conj_t conjx, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* x, inc_t incx, \
             ctype* a, inc_t rs_a, inc_t cs_a  \
     ) \
{ \
    bli_##ch##syr( uploa, conjx, m, alpha, x, incx, a, rs_a, cs_a ); \
} \
\
inline void bli_her2 \
     ( \
             uplo_t uploa, \
             conj_t conjx, \
             conj_t conjy, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* x, inc_t incx, \
       const ctype* y, inc_t incy, \
             ctype* a, inc_t rs_a, inc_t cs_a  \
     ) \
{ \
    bli_##ch##her2( uploa, conjx, conjy, m, alpha, x, incx, y, incy, a, rs_a, cs_a ); \
} \
\
inline void bli_syr2 \
     ( \
             uplo_t uploa, \
             conj_t conjx, \
             conj_t conjy, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* x, inc_t incx, \
       const ctype* y, inc_t incy, \
             ctype* a, inc_t rs_a, inc_t cs_a  \
     ) \
{ \
    bli_##ch##syr2( uploa, conjx, conjy, m, alpha, x, incx, y, incy, a, rs_a, cs_a ); \
} \
\
inline void bli_shr2 \
     ( \
             uplo_t uploa, \
             conj_t conjx, \
             conj_t conjy, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* x, inc_t incx, \
       const ctype* y, inc_t incy, \
             ctype* a, inc_t rs_a, inc_t cs_a  \
     ) \
{ \
    bli_##ch##shr2( uploa, conjx, conjy, m, alpha, x, incx, y, incy, a, rs_a, cs_a ); \
} \
\
inline void bli_skr2 \
     ( \
             uplo_t uploa, \
             conj_t conjx, \
             conj_t conjy, \
             dim_t  m, \
       const ctype* alpha, \
       const ctype* x, inc_t incx, \
       const ctype* y, inc_t incy, \
             ctype* a, inc_t rs_a, inc_t cs_a  \
     ) \
{ \
    bli_##ch##skr2( uploa, conjx, conjy, m, alpha, x, incx, y, incy, a, rs_a, cs_a ); \
} \
\
inline void bli_trmv \
     ( \
             uplo_t  uploa, \
             trans_t transa, \
             diag_t  diaga, \
             dim_t   m, \
       const ctype*  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
             ctype*  x, inc_t incx  \
     ) \
{ \
    bli_##ch##trmv( uploa, transa, diaga, m, alpha, a, rs_a, cs_a, x, incx ); \
} \
\
inline void bli_trsv \
     ( \
             uplo_t  uploa, \
             trans_t transa, \
             diag_t  diaga, \
             dim_t   m, \
       const ctype*  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
             ctype*  x, inc_t incx  \
     ) \
{ \
    bli_##ch##trsv( uploa, transa, diaga, m, alpha, a, rs_a, cs_a, x, incx ); \
}

MARRAY_FOR_EACH_TYPE

/******************************************************************************
 *
 * Level 3 BLAS, C++ overloads
 *
 *****************************************************************************/
#undef MARRAY_FOR_EACH_TYPE_BODY
#define MARRAY_FOR_EACH_TYPE_BODY(ctype, ctypef, ctyper, CH, ch, CHR, chr) \
\
inline void bli_gemm \
     ( \
             trans_t transa, \
             trans_t transb, \
             dim_t   m, \
             dim_t   n, \
             dim_t   k, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##gemm( transa, transb, m, n, k, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_hemm \
     ( \
             side_t  side, \
             uplo_t  uploa, \
             conj_t  conja, \
             trans_t transb, \
             dim_t   m, \
             dim_t   n, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##hemm( side, uploa, conja, transb, m, n, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_symm \
     ( \
             side_t  side, \
             uplo_t  uploa, \
             conj_t  conja, \
             trans_t transb, \
             dim_t   m, \
             dim_t   n, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##symm( side, uploa, conja, transb, m, n, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_shmm \
     ( \
             side_t  side, \
             uplo_t  uploa, \
             conj_t  conja, \
             trans_t transb, \
             dim_t   m, \
             dim_t   n, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##shmm( side, uploa, conja, transb, m, n, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_skmm \
     ( \
             side_t  side, \
             uplo_t  uploa, \
             conj_t  conja, \
             trans_t transb, \
             dim_t   m, \
             dim_t   n, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##skmm( side, uploa, conja, transb, m, n, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_herk \
     ( \
             uplo_t   uploc, \
             trans_t  transa, \
             dim_t    m, \
             dim_t    k, \
       const ctyper&  alpha, \
       const ctype*   a, inc_t rs_a, inc_t cs_a, \
       const ctyper&  beta, \
             ctype*   c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##herk( uploc, transa, m, k, &alpha, a, rs_a, cs_a, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_her2k \
     ( \
             uplo_t   uploc, \
             trans_t  transa, \
             trans_t  transb, \
             dim_t    m, \
             dim_t    k, \
       const ctype&   alpha, \
       const ctype*   a, inc_t rs_a, inc_t cs_a, \
       const ctype*   b, inc_t rs_b, inc_t cs_b, \
       const ctyper&  beta, \
             ctype*   c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##her2k( uploc, transa, transb, m, k, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_syrk \
     ( \
             uplo_t  uploc, \
             trans_t transa, \
             dim_t   m, \
             dim_t   k, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##syrk( uploc, transa, m, k, &alpha, a, rs_a, cs_a, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_syr2k \
     ( \
             uplo_t  uploc, \
             trans_t transa, \
             trans_t transb, \
             dim_t   m, \
             dim_t   k, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##syr2k( uploc, transa, transb, m, k, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_shr2k \
     ( \
             uplo_t   uploc, \
             trans_t  transa, \
             trans_t  transb, \
             dim_t    m, \
             dim_t    k, \
       const ctype&   alpha, \
       const ctype*   a, inc_t rs_a, inc_t cs_a, \
       const ctype*   b, inc_t rs_b, inc_t cs_b, \
       const ctyper&  beta, \
             ctype*   c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##shr2k( uploc, transa, transb, m, k, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_skr2k \
     ( \
             uplo_t  uploc, \
             trans_t transa, \
             trans_t transb, \
             dim_t   m, \
             dim_t   k, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##skr2k( uploc, transa, transb, m, k, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_gemmt \
     ( \
             uplo_t  uploc, \
             trans_t transa, \
             trans_t transb, \
             dim_t   m, \
             dim_t   k, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##gemmt( uploc, transa, transb, m, k, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_trmm3 \
     ( \
             side_t  side, \
             uplo_t  uploa, \
             trans_t transa, \
             diag_t  diaga, \
             trans_t transb, \
             dim_t   m, \
             dim_t   n, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
       const ctype*  b, inc_t rs_b, inc_t cs_b, \
       const ctype&  beta, \
             ctype*  c, inc_t rs_c, inc_t cs_c  \
     ) \
{ \
    bli_##ch##trmm3( side, uploa, transa, diaga, transb, m, n, &alpha, a, rs_a, cs_a, b, rs_b, cs_b, &beta, c, rs_c, cs_c ); \
} \
\
inline void bli_trmm \
     ( \
             side_t  side, \
             uplo_t  uploa, \
             trans_t transa, \
             diag_t  diaga, \
             dim_t   m, \
             dim_t   n, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
             ctype*  b, inc_t rs_b, inc_t cs_b  \
     ) \
{ \
    bli_##ch##trmm( side, uploa, transa, diaga, m, n, &alpha, a, rs_a, cs_a, b, rs_b, cs_b ); \
} \
\
inline void bli_trsm \
     ( \
             side_t  side, \
             uplo_t  uploa, \
             trans_t transa, \
             diag_t  diaga, \
             dim_t   m, \
             dim_t   n, \
       const ctype&  alpha, \
       const ctype*  a, inc_t rs_a, inc_t cs_a, \
             ctype*  b, inc_t rs_b, inc_t cs_b  \
     ) \
{ \
    bli_##ch##trsm( side, uploa, transa, diaga, m, n, &alpha, a, rs_a, cs_a, b, rs_b, cs_b ); \
}

MARRAY_FOR_EACH_TYPE

#endif //BLIS_H

/******************************************************************************
 *
 * Level 1 BLAS, MArray wrappers
 *
 *****************************************************************************/
template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1>>
swapv(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());

#ifdef BLIS_H
    bli_swapv(x.length(), x.data(), x.stride(), y.data(), y.stride());
#else
    swap(x.length(), x.data(), x.stride(), y.data(), y.stride());
#endif
}

template <typename T>
std::enable_if_t<detail::is_marray_like_v<T,1>>
conj(T&& x_)
{
    auto x = x_.view();

#ifdef BLIS_H
    bli_scalv(BLIS_CONJUGATE, x.length(), detail::value_type<T>{1}, x.data(), x.stride());
#else
    auto n = x.length();
    auto ptr = x.data();
    auto stride = x.stride();
    for (len_type i = 0;i < n;i++)
        ptr[i*stride] = std::conj(ptr[i*stride]);
#endif
}

template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<U,1>>
scal(T alpha, U&& x_)
{
    auto x = x_.view();

#ifdef BLIS_H
    bli_scalv(BLIS_NO_CONJUGATE, x.length(), alpha, x.data(), x.stride());
#else
    scal(x.length(), alpha, x.data(), x.stride());
#endif
}

template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1>>
copy(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());

#ifdef BLIS_H
    bli_copyv(BLIS_NO_CONJUGATE, x.length(), x.data(), x.stride(), y.data(), y.stride());
#else
    copy(x.length(), x.data(), x.stride(), y.data(), y.stride());
#endif
}

template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1>>
axpy(T alpha, U&& x_, V&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());

#ifdef BLIS_H
    bli_axpyv(BLIS_NO_CONJUGATE, x.length(), alpha, x.data(), x.stride(), y.data(), y.stride());
#else
    axpy(x.length(), alpha, x.data(), x.stride(), y.data(), y.stride());
#endif
}

template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto dotu(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());

#ifdef BLIS_H
    return bli_dotv(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, x.length(), x.data(), x.stride(), y.data(), y.stride());
#else
    return dotu(x.length(), x.data(), x.stride(), y.data(), y.stride());
#endif
}

template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto dotc(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());

#ifdef BLIS_H
    return bli_dotv(BLIS_NO_CONJUGATE, BLIS_CONJUGATE, x.length(), x.data(), x.stride(), y.data(), y.stride());
#else
    return dotc(x.length(), x.data(), x.stride(), y.data(), y.stride());
#endif
}

template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto dot(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());

#ifdef BLIS_H
    return bli_dotv(BLIS_NO_CONJUGATE, BLIS_CONJUGATE, x.length(), x.data(), x.stride(), y.data(), y.stride());
#else
    return dot(x.length(), x.data(), x.stride(), y.data(), y.stride());
#endif
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto nrm2(T&& x_)
{
    auto x = x_.view();

#ifdef BLIS_H
    return bli_normfv(x.length(), x.data(), x.stride());
#else
    return nrm2(x.length(), x.data(), x.stride());
#endif
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto asum(T&& x_)
{
    auto x = x_.view();

#ifdef BLIS_H
    return bli_norm1v(x.length(), x.data(), x.stride());
#else
    return asum(x.length(), x.data(), x.stride());
#endif
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto iamax(T&& x_)
{
    auto x = x_.view();

#ifdef BLIS_H
    return bli_amaxv(x.length(), x.data(), x.stride());
#else
    return iamax(x.length(), x.data(), x.stride());
#endif
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto amax(T&& x_)
{
    auto x = x_.view();

#ifdef BLIS_H
    return bli_normiv(x.length(), x.data(), x.stride());
#else
    return amax(x.length(), x.data(), x.stride());
#endif
}

/******************************************************************************
 *
 * Level 2 BLAS, MArray wrappers
 *
 *****************************************************************************/

/**
 * Perform the matrix-vector multiplication \f$ y = \alpha Ax + \beta y \f$.
 *
 * The "left-side" multiplication \f$ y = \alpha xA + \beta y \f$ can be performed using
 * `gemv(alpha, A.T(), x, beta, y)`.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `n`.
 *
 * @param beta  Scalar factor for the original vector `y`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,1>>
gemv(T alpha, U&& A, V&& x, W beta, X&& y)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);
    auto n = A_.length(1);

    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(x_.length() == n);

#ifdef BLIS_H
    bli_gemv(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, m, n,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    x_.data(), x_.stride(),
              beta, y_.data(), y_.stride());
#else
    char transa = 'N';

    if (A_.stride(0) > 1)
    {
        A_.transpose();
        transa = 'T';
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    gemv(transa, A_.length(0), A_.length(1),
         alpha, A_.data(), A_.stride(1),
                x_.data(), x_.stride(),
          beta, y_.data(), y_.stride());
#endif
}

/**
 * Perform the matrix-vector multiplication \f$ y = Ax \f$.
 *
 * The "left-side" multiplication \f$ y = xA \f$ can be performed using
 * `gemv(A.T(), x, y)`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `n`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1>>
gemv(T&& A, U&& x, V&& y)
{
    gemv(1.0, A, x, 0.0, y);
}

/**
 * Return the result of the matrix-vector multiplication \f$ y = \alpha Ax \f$.
 *
 * The "left-side" multiplication \f$ y = \alpha xA \f$ can be performed using
 * `gemv(alpha, A.T(), x)`.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `n`.
 *
 * @return      A vector of length `m` holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,1>>>
auto gemv(T alpha, U&& A, V&& x)
{
    marray<detail::value_type<U>,1> y({A.length(0)}, uninitialized);
    gemv(alpha, A, x, 0.0, y);
    return y;
}

/**
 * Return the result of the matrix-vector multiplication \f$ y = Ax \f$.
 *
 * The "left-side" multiplication \f$ y = xA \f$ can be performed using
 * `gemv(A.T(), x)`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `n`.
 *
 * @return      A vector of length `m` holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,1>>>
auto gemv(T&& A, U&& x)
{
    return gemv(1.0, A, x);
}

/**
 * Perform the Hermitian matrix-vector multiplication \f$ y = \alpha Ax + \beta y \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original vector `y`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,1>>
hemv(char uplo, T alpha, U&& A, V&& x, W beta, X&& y)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    bli_hemv(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    x_.data(), x_.stride(),
              beta, y_.data(), y_.stride());
#else
    if (A_.stride(0) > 1)
    {
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;

        if constexpr (detail::is_complex_v<detail::value_type<V>>)
        {
            // Compute (alpha Ax + beta y)^H = alpha' x^H A + beta' y^H

            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<detail::value_type<V>,1> conjx(x_);

            alpha = std::conj(alpha);
            beta = std::conj(beta);
            conj(conjx);
            conj(y_);

            hemv(uplo, m,
                 alpha, A_.data(), A_.stride(1),
                        conjx.data(), conjx.stride(),
                  beta, y_.data(), y_.stride());

            conj(y_);

            return;
        }
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    hemv(uplo, m,
         alpha, A_.data(), A_.stride(1),
                x_.data(), x_.stride(),
          beta, y_.data(), y_.stride());
#endif
}

/**
 * Perform the Hermitian matrix-vector multiplication \f$ y = Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1>>
hemv(char uplo, T&& A, U&& x, V&& y)
{
    hemv(uplo, 1.0, A, x, 0.0, y);
}

/**
 * Return the result of the Hermitian matrix-vector multiplication \f$ y = \alpha Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,1>>>
auto hemv(char uplo, T alpha, U&& A, V&& x)
{
    marray<detail::value_type<U>,1> y({A.length(0)}, uninitialized);
    hemv(uplo, alpha, A, x, 0.0, y);
    return y;
}

/**
 * Return the result of the Hermitian matrix-vector multiplication \f$ y = Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,1>>>
auto hemv(char uplo, T&& A, U&& x)
{
    return hemv(uplo, 1.0, A, x);
}

/**
 * Perform the triangular matrix-vector multiplication \f$ y = \alpha Ax + \beta y \f$.
 *
 * The "left-side" multiplication \f$ y = \alpha xA + \beta y \f$ can be performed using
 * `trmv(fliplu(uplo), diag, alpha, A.T(), x, beta, y)`, where fliplu('L') == 'U' and
 * vice versa.
 *
 * @param uplo  'L' if A is lower-triangular or 'U' if A is upper-triangular.
 *
 * @param diag  'U' if A is unit-diagonal, or 'N' if A is non-unit-diagonal.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` triangular matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original vector `y`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,1>>
trmv(char uplo, char diag, T alpha, U&& A, V&& x, W beta, X&& y)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

    marray<detail::value_type<X>,1> y0;

    if (beta != 0.0)
        y0.reset(y);

    copy(x, y);

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    MARRAY_ASSERT(diag == 'U' || diag == 'N');
    bli_trmv(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_TRANSPOSE,
             diag == 'U' ? BLIS_UNIT_DIAG : BLIS_NONUNIT_DIAG, m,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    y_.data(), y_.stride());
#else
    char transa = 'N';

    if (A_.stride(0) > 1)
    {
        transa = 'T';
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    trmv(uplo, transa, diag, m,
         alpha, A_.data(), A_.stride(1),
                y_.data(), y_.stride());
#endif

    if (beta != 0.0)
        axpy(beta, y0, y);
}

/**
 * Perform the triangular matrix-vector multiplication \f$ y = Ax \f$.
 *
 * The "left-side" multiplication \f$ y = xA \f$ can be performed using
 * `trmv(fliplu(uplo), diag, A.T(), x, y)`, where fliplu('L') == 'U' and
 * vice versa.
 *
 * @param uplo  'L' if A is lower-triangular or 'U' if A is upper-triangular.
 *
 * @param diag  'U' if A is unit-diagonal, or 'N' if A is non-unit-diagonal.
 *
 * @param A     A `m`x`m` triangular matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1>>
trmv(char uplo, char diag, T&& A, U&& x, V&& y)
{
    trmv(uplo, diag, 1.0, A, x, 0.0, y);
}

/**
 * Return the result of the triangular matrix-vector multiplication \f$ y = \alpha Ax \f$.
 *
 * The "left-side" multiplication \f$ \alpha xA \f$ can be performed using
 * `trmv(fliplu(uplo), diag, alpha, A.T(), x)`, where fliplu('L') == 'U' and
 * vice versa.
 *
 * @param uplo  'L' if A is lower-triangular or 'U' if A is upper-triangular.
 *
 * @param diag  'U' if A is unit-diagonal, or 'N' if A is non-unit-diagonal.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` triangular matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,1>>>
auto trmv(char uplo, char diag, T alpha, U&& A, V&& x)
{
    marray<detail::value_type<U>,1> y({A.length(0)}, uninitialized);
    trmv(uplo, diag, alpha, A, x, 0.0, y);
    return y;
}

/**
 * Return the result of the triangular matrix-vector multiplication \f$ y = Ax \f$.
 *
 * The "left-side" multiplication \f$ xA \f$ can be performed using
 * `trmv(fliplu(uplo), diag, A.T(), x)`, where fliplu('L') == 'U' and
 * vice versa.
 *
 * @param uplo  'L' if A is lower-triangular or 'U' if A is upper-triangular.
 *
 * @param diag  'U' if A is unit-diagonal, or 'N' if A is non-unit-diagonal.
 *
 * @param A     A `m`x`m` triangular matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,1>>>
auto trmv(char uplo, char diag, T&& A, U&& x)
{
    return trmv(uplo, diag, 1.0, A, x);
}

/**
 * Perform a triangular solve \f$ b = Ax \f$, where b is originally stored in the solution vector x.
 *
 * The "left-side" solve \f$ b = xA \f$ can be performed using
 * `trsv(fliplu(uplo), diag, A.T(), x)`, where fliplu('L') == 'U' and
 * vice versa.
 *
 * @param uplo  'L' if A is lower-triangular or 'U' if A is upper-triangular.
 *
 * @param diag  'U' if A is unit-diagonal, or 'N' if A is non-unit-diagonal.
 *
 * @param A     A `m`x`m` triangular matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,1>>
trsv(char uplo, char diag, T&& A, U&& x)
{
    auto A_ = A.view();
    auto x_ = x.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    MARRAY_ASSERT(diag == 'U' || diag == 'N');
    bli_trsv(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_TRANSPOSE,
             diag == 'U' ? BLIS_UNIT_DIAG : BLIS_NONUNIT_DIAG, m,
             A_.data(), A_.stride(0), A_.stride(1),
             x_.data(), x_.stride());
#else
    char transa = 'N';

    if (A_.stride(0) > 1)
    {
        transa = 'T';
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    trsv(uplo, transa, diag, m,
         A_.data(), A_.stride(1),
         x_.data(), x_.stride());
#endif
}

/**
 * Perform the Hermitian outer product \f$ A = \alpha xx^\text{H} + \beta A \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xx^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original matrix `A`.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<W,2>>
__attribute__((always_inline))
her(char uplo, T alpha, U&& x, V beta, W&& A)
{
    auto A_ = A.view();
    auto x_ = x.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

    if (beta == 0.0) A_ = 0;
    else if (beta != 1.0) A_ *= beta;

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    bli_her(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_CONJUGATE, m,
            alpha, x_.data(), x_.stride(),
                   A_.data(), A_.stride(0), A_.stride(1));
#else
    if (A_.stride(0) > 1 && 0)
    {
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;

        if constexpr (detail::is_complex_v<detail::value_type<U>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<detail::value_type<U>,1> conjx(x_);
            conj(conjx);

            her(uplo, m,
                alpha, conjx.data(), conjx.stride(),
                       A_.data(), A_.stride(1));

            return;
        }
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    her(uplo, m,
        alpha, x_.data(), x_.stride(),
               A_.data(), A_.stride(1));
#endif
}

/**
 * Perform the Hermitian outer product \f$ A = xx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,2>>
her(char uplo, T&& x, U&& A)
{
    her(uplo, 1.0, x, 0.0, A);
}

/**
 * Return the result of the Hermitian outer product \f$ A = \alpha xx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xx^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<U,1>>>
auto her(char uplo, T alpha, U&& x)
{
    marray<detail::value_type<U>,2> A({x.length(), x.length()}, uninitialized);
    her(uplo, alpha, x, 0.0, A);
    return A;
}

/**
 * Return the result of the Hermitian outer product \f$ A = xx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto her(char uplo, T&& x)
{
    return her(uplo, 1.0, x);
}

/**
 * Perform the Hermitian outer product \f$ A = \alpha xy^\text{H} + \alpha' yx^\text{H} + \beta A \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original matrix `A`.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,2>>
her2(char uplo, T alpha, U&& x, V&& y, W beta, X&& A)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

    if (beta == 0.0) A_ = 0;
    else if (beta != 1.0) A_ *= beta;

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    bli_her2(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m,
             alpha, x_.data(), x_.stride(),
                    y_.data(), y_.stride(),
                    A_.data(), A_.stride(0), A_.stride(1));
#else
    if (A_.stride(0) > 1)
    {
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;

        if constexpr (detail::is_complex_v<detail::value_type<U>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<detail::value_type<U>,1> conjx(x_);
            marray<detail::value_type<V>,1> conjy(y_);
            conj(conjx);
            conj(conjy);

            her2(uplo, m,
                 alpha, conjx.data(), conjx.stride(),
                        conjy.data(), conjy.stride(),
                        A_.data(), A_.stride(1));

            return;
        }
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    her2(uplo, m,
         alpha, x_.data(), x_.stride(),
                y_.data(), y_.stride(),
                A_.data(), A_.stride(1));
#endif
}

/**
 * Perform the Hermitian outer product \f$ A = xy^\text{H} + yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 *
 * @param A     A `m`x`m` Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,2>>
her2(char uplo, T&& x, U&& y, V&& A)
{
    her2(uplo, 1.0, x, y, 0.0, A);
}

/**
 * Return the result of the Hermitian outer product \f$ A = \alpha xy^\text{H} + \alpha' yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,1> &&
                     detail::is_marray_like_v<V,1>>>
auto her2(char uplo, T alpha, U&& x, V&& y)
{
    marray<detail::value_type<U>,2> A({x.length(), x.length()}, uninitialized);
    her2(uplo, alpha, x, y, 0.0, A);
    return A;
}

/**
 * Return the result of the Hermitian outer product \f$ A = xy^\text{H} + yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto her2(char uplo, T&& x, U&& y)
{
    return her2(uplo, 1.0, x, y);
}

/**
 * Perform the outer product \f$ A = \alpha xy^\text{H} + \beta A \f$.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @param beta  Scalar factor for the original matrix `A`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,2>>
ger(T alpha, U&& x, V&& y, W beta, X&& A)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);
    auto n = A_.length(1);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(y_.length() == n);

    if (beta == 0.0) A_ = 0;
    else if (beta != 1.0) A_ *= beta;

#ifdef BLIS_H
    bli_ger(BLIS_NO_CONJUGATE, BLIS_CONJUGATE, m, n,
            alpha, x_.data(), x_.stride(),
                   y_.data(), y_.stride(),
                   A_.data(), A_.stride(0), A_.stride(1));
#else
    if (A_.stride(0) > 1)
    {
        A_.transpose();
        x_.swap(y_);

        if constexpr (detail::is_complex_v<detail::value_type<V>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<detail::value_type<V>,1> conjx(x_);
            conj(conjx);

            geru(A_.length(0), A_.length(1),
                 alpha, x_.data(), x_.stride(),
                        y_.data(), y_.stride(),
                        A_.data(), A_.stride(1));

            return;
        }
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    gerc(A_.length(0), A_.length(1),
         alpha, x_.data(), x_.stride(),
                y_.data(), y_.stride(),
                A_.data(), A_.stride(1));
#endif
}

/**
 * Perform the outer product \f$ A = xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,2>>
ger(T&& x, U&& y, V&& A)
{
    ger(1.0, x, y, 0.0, A);
}

/**
 * Return the result of the outer product \f$ A = \alpha xy^\text{H} \f$.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @return      A `m`x`n` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,1> &&
                     detail::is_marray_like_v<V,1>>>
auto ger(T alpha, U&& x, V&& y)
{
    marray<detail::value_type<U>,2> A({x.length(), y.length()}, uninitialized);
    ger(alpha, x, y, 0.0, A);
    return A;
}

/**
 * Return the result of the outer product \f$ A = xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @return      A `m`x`n` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto ger(T&& x, U&& y)
{
    return ger(1.0, x, y);
}

/**
 * Perform the outer product \f$ A = \alpha xy^\text{H} + \beta A \f$.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @param beta  Scalar factor for the original matrix `A`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,2>>
gerc(T alpha, U&& x, V&& y, W beta, X&& A)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);
    auto n = A_.length(1);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(y_.length() == n);

    if (beta == 0.0) A_ = 0;
    else if (beta != 1.0) A_ *= beta;

#ifdef BLIS_H
    bli_ger(BLIS_NO_CONJUGATE, BLIS_CONJUGATE, m, n,
            alpha, x_.data(), x_.stride(),
                   y_.data(), y_.stride(),
                   A_.data(), A_.stride(0), A_.stride(1));
#else
    if (A_.stride(0) > 1)
    {
        A_.transpose();
        x_.swap(y_);

        if constexpr (detail::is_complex_v<detail::value_type<V>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<detail::value_type<V>,1> conjx(x_);
            conj(conjx);

            geru(A_.length(0), A_.length(1),
                 alpha, x_.data(), x_.stride(),
                        y_.data(), y_.stride(),
                        A_.data(), A_.stride(1));

            return;
        }
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    gerc(A_.length(0), A_.length(1),
         alpha, x_.data(), x_.stride(),
                y_.data(), y_.stride(),
                A_.data(), A_.stride(1));
#endif
}

/**
 * Perform the outer product \f$ A = xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,2>>
gerc(T&& x, U&& y, V&& A)
{
    gerc(1.0, x, y, 0.0, A);
}

/**
 * Return the result of the outer product \f$ A = \alpha xy^\text{H} \f$.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @return      A `m`x`n` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,1> &&
                     detail::is_marray_like_v<V,1>>>
auto gerc(T alpha, U&& x, V&& y)
{
    marray<detail::value_type<U>,2> A({x.length(), y.length()}, uninitialized);
    gerc(alpha, x, y, 0.0, A);
    return A;
}

/**
 * Return the result of the outer product \f$ A = xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @return      A `m`x`n` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto gerc(T&& x, U&& y)
{
    return gerc(1.0, x, y);
}

/**
 * Perform the outer product \f$ A = \alpha xy^\text{T} + \beta A \f$.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{T} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @param beta  Scalar factor for the original matrix `A`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,2>>
geru(T alpha, U&& x, V&& y, W beta, X&& A)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);
    auto n = A_.length(1);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(y_.length() == n);

    if (beta == 0.0) A_ = 0;
    else if (beta != 1.0) A_ *= beta;

#ifdef BLIS_H
    bli_ger(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m, n,
            alpha, x_.data(), x_.stride(),
                   y_.data(), y_.stride(),
                   A_.data(), A_.stride(0), A_.stride(1));
#else
    if (A_.stride(0) > 1)
    {
        A_.transpose();
        x_.swap(y_);
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    geru(A_.length(0), A_.length(1),
         alpha, x_.data(), x_.stride(),
                y_.data(), y_.stride(),
                A_.data(), A_.stride(1));
#endif
}

/**
 * Perform the outer product \f$ A = xy^\text{T} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @param A     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,2>>
geru(T&& x, U&& y, V&& A)
{
    geru(1.0, x, y, 0.0, A);
}

/**
 * Return the result of the outer product \f$ A = \alpha xy^\text{T} \f$.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{T} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @return      A `m`x`n` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,1> &&
                     detail::is_marray_like_v<V,1>>>
auto geru(T alpha, U&& x, V&& y)
{
    marray<detail::value_type<U>,2> A({x.length(), y.length()}, uninitialized);
    geru(alpha, x, y, 0.0, A);
    return A;
}

/**
 * Return the result of the outer product \f$ A = xy^\text{T} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `n`.
 *
 * @return      A `m`x`n` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto geru(T&& x, U&& y)
{
    return geru(1.0, x, y);
}

#ifdef BLIS_H

/**
 * Perform the skew-Hermitian matrix-vector multiplication \f$ y = \alpha Ax + \beta y \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` skew-Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original vector `y`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,1>>
shmv(char uplo, T alpha, U&& A, V&& x, W beta, X&& y)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    bli_shmv(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    x_.data(), x_.stride(),
              beta, y_.data(), y_.stride());
}

/**
 * Perform the skew-Hermitian matrix-vector multiplication \f$ y = Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` skew-Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1>>
shmv(char uplo, T&& A, U&& x, V&& y)
{
    shmv(uplo, 1.0, A, x, 0.0, y);
}

/**
 * Return the result of the skew-Hermitian matrix-vector multiplication \f$ y = \alpha Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` skew-Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,1>>>
auto shmv(char uplo, T alpha, U&& A, V&& x)
{
    marray<detail::value_type<U>,1> y({A.length(0)}, uninitialized);
    shmv(uplo, alpha, A, x, 0.0, y);
    return y;
}

/**
 * Return the result of the skew-Hermitian matrix-vector multiplication \f$ y = Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` skew-Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,1>>>
auto shmv(char uplo, T&& A, U&& x)
{
    return shmv(uplo, 1.0, A, x);
}

/**
 * Perform the skew-symmetric matrix-vector multiplication \f$ y = \alpha Ax + \beta y \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` skew-symmetric matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original vector `y`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,1>>
skmv(char uplo, T alpha, U&& A, V&& x, W beta, X&& y)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    bli_skmv(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    x_.data(), x_.stride(),
              beta, y_.data(), y_.stride());
}

/**
 * Perform the skew-symmetric matrix-vector multiplication \f$ y = Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` skew-symmetric matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1>>
skmv(char uplo, T&& A, U&& x, V&& y)
{
    skmv(uplo, 1.0, A, x, 0.0, y);
}

/**
 * Return the result of the skew-symmetric matrix-vector multiplication \f$ y = \alpha Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `Ax`.
 *
 * @param A     A `m`x`m` skew-symmetric matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,1>>>
auto skmv(char uplo, T alpha, U&& A, V&& x)
{
    marray<detail::value_type<U>,1> y({A.length(0)}, uninitialized);
    skmv(uplo, alpha, A, x, 0.0, y);
    return y;
}

/**
 * Return the result of the skew-symmetric matrix-vector multiplication \f$ y = Ax \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` skew-symmetric matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param x     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,1>>>
auto skmv(char uplo, T&& A, U&& x)
{
    return skmv(uplo, 1.0, A, x);
}

/**
 * Perform the skew-Hermitian outer product \f$ A = \alpha xy^\text{H} - \alpha' yx^\text{H} + \beta A \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original matrix `A`.
 *
 * @param A     A `m`x`m` skew-Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,2>>
shr2(char uplo, T alpha, U&& x, V&& y, W beta, X&& A)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

    if (beta == 0.0) A_ = 0;
    else if (beta != 1.0) A_ *= beta;

    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    bli_shr2(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m,
             alpha, x_.data(), x_.stride(),
                    y_.data(), y_.stride(),
                    A_.data(), A_.stride(0), A_.stride(1));
}

/**
 * Perform the skew-Hermitian outer product \f$ A = xy^\text{H} - yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 *
 * @param A     A `m`x`m` skew-Hermitian matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,2>>
shr2(char uplo, T&& x, U&& y, V&& A)
{
    shr2(uplo, 1.0, x, y, 0.0, A);
}

/**
 * Return the result of the skew-Hermitian outer product \f$ A = \alpha xy^\text{H} - \alpha' yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,1> &&
                     detail::is_marray_like_v<V,1>>>
auto shr2(char uplo, T alpha, U&& x, V&& y)
{
    marray<detail::value_type<U>,2> A({x.length(), x.length()}, uninitialized);
    shr2(uplo, alpha, x, y, 0.0, A);
    return A;
}

/**
 * Return the result of the skew-Hermitian outer product \f$ A = xy^\text{H} - yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto shr2(char uplo, T&& x, U&& y)
{
    return shr2(uplo, 1.0, x, y);
}

/**
 * Perform the skew-symmetric outer product \f$ A = \alpha xy^\text{H} - \alpha' yx^\text{H} + \beta A \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 *
 * @param beta  Scalar factor for the original matrix `A`.
 *
 * @param A     A `m`x`m` skew-symmetric matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1> &&
                 detail::is_marray_like_v<X,2>>
skr2(char uplo, T alpha, U&& x, V&& y, W beta, X&& A)
{
    auto A_ = A.view();
    auto x_ = x.view();
    auto y_ = y.view();

    auto m = A_.length(0);

    MARRAY_ASSERT(x_.length() == m);
    MARRAY_ASSERT(y_.length() == m);
    MARRAY_ASSERT(A_.length(1) == m);

    if (beta == 0.0) A_ = 0;
    else if (beta != 1.0) A_ *= beta;

    MARRAY_ASSERT(uplo == 'U' || uplo == 'L');
    bli_skr2(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m,
             alpha, x_.data(), x_.stride(),
                    y_.data(), y_.stride(),
                    A_.data(), A_.stride(0), A_.stride(1));
}

/**
 * Perform the skew-symmetric outer product \f$ A = xy^\text{H} - yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 *
 * @param A     A `m`x`m` skew-symmetric matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,2>>
skr2(char uplo, T&& x, U&& y, V&& A)
{
    skr2(uplo, 1.0, x, y, 0.0, A);
}

/**
 * Return the result of the skew-symmetric outer product \f$ A = \alpha xy^\text{H} - \alpha' yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product \f$ xy^\text{H} \f$.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,1> &&
                     detail::is_marray_like_v<V,1>>>
auto skr2(char uplo, T alpha, U&& x, V&& y)
{
    marray<detail::value_type<U>,2> A({x.length(), x.length()}, uninitialized);
    skr2(uplo, alpha, x, y, 0.0, A);
    return A;
}

/**
 * Return the result of the skew-symmetric outer product \f$ A = xy^\text{H} - yx^\text{H} \f$.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param x     A vector or vector view of length `m`.
 *
 * @param y     A vector or vector view of length `m`.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto skr2(char uplo, T&& x, U&& y)
{
    return skr2(uplo, 1.0, x, y);
}

#endif

/******************************************************************************
 *
 * Level 3 BLAS, MArray wrappers
 *
 *****************************************************************************/

/**
 * Perform the matrix multiplication \f$ C = \alpha AB + \beta C \f$.
 *
 * @param alpha Scalar factor for the product `AB`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `k`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<X,2>>
gemm(T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = C_.length(0);
    auto n = C_.length(1);
    auto k = A_.length(1);

    MARRAY_ASSERT(A_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == n);
    MARRAY_ASSERT(B_.length(0) == k);

#ifdef BLIS_H
    bli_gemm(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, m, n, k,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    B_.data(), B_.stride(0), B_.stride(1),
              beta, C_.data(), C_.stride(0), C_.stride(1));
#else
    if (C.stride(0) > 1)
    {
        std::swap(m, n);
        A_.transpose();
        B_.transpose();
        C_.transpose();
        A_.swap(B_);
    }

    char transa = A_.stride(0) > 1 ? 'T' : 'N';
    char transb = B_.stride(0) > 1 ? 'T' : 'N';

    if (transa == 'T') A_.transpose();
    if (transb == 'T') B_.transpose();

    MARRAY_ASSERT(A_.stride(0) == 1 && A_.stride(1) > 1);
    MARRAY_ASSERT(B_.stride(0) == 1 && B_.stride(1) > 1);
    MARRAY_ASSERT(C_.stride(0) == 1 && C_.stride(1) > 1);

    gemm(transa, transb, m, n, k,
         alpha, A_.data(), A_.stride(1),
                B_.data(), B_.stride(1),
          beta, C_.data(), C_.stride(1));
#endif
}

/**
 * Perform the matrix multiplication \f$ C = AB \f$.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `k`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
gemm(T&& A, U&& B, V&& C)
{
    gemm(1.0, A, B, 0.0, C);
}

/**
 * Return the result of the matrix multiplication \f$ C = \alpha AB \f$.
 *
 * @param alpha Scalar factor for the product `AB`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `k`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`n` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto gemm(T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({A.length(0), B.length(1)}, uninitialized);
    gemm(alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the matrix multiplication \f$ C = AB \f$.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `k`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`n` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto gemm(T&& A, U&& B)
{
    return gemm(1.0, A, B);
}

/**
 * Perform the matrix multiplication \f$ C = \alpha AB + \beta C \f$ or
 * \f$ C = \alpha BA + \beta C \f$, where A is a symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
symm(char side, char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = C_.length(0);
    auto n = C_.length(1);

    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == n);
    MARRAY_ASSERT(A_.length(0) == side == 'L' ? m : n);
    MARRAY_ASSERT(A_.length(1) == side == 'L' ? m : n);

#ifdef BLIS_H
    MARRAY_ASSERT(side == 'L' || side == 'R');
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_symm(side == 'L' ? BLIS_LEFT : BLIS_RIGHT, uplo == 'U' ? BLIS_UPPER : BLIS_LOWER,
             BLIS_NO_CONJUGATE, BLIS_NO_TRANSPOSE, m, n,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    B_.data(), B_.stride(0), B_.stride(1),
              beta, C_.data(), C_.stride(0), C_.stride(1));
#else
    if (C.stride(0) > 1)
    {
        side = side == 'L' ? 'R' :
               side == 'R' ? 'L' : side;
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        std::swap(m, n);
        A_.transpose();
        B_.transpose();
        C_.transpose();
    }

    MARRAY_ASSERT(A_.stride(0) == 1 && A_.stride(1) > 1);
    MARRAY_ASSERT(B_.stride(0) == 1 && B_.stride(1) > 1);
    MARRAY_ASSERT(C_.stride(0) == 1 && C_.stride(1) > 1);

    symm(side, uplo, m, n,
         alpha, A_.data(), A_.stride(1),
                B_.data(), B_.stride(1),
          beta, C_.data(), C_.stride(1));
#endif
}

/**
 * Perform the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
symm(char side, char uplo, T&& A, U&& B, V&& C)
{
    symm(side, uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the matrix multiplication \f$ C = \alpha AB \f$ or
 * \f$ C = \alpha BA \f$, where A is a symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto symm(char side, char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({B.length(0), B.length(1)}, uninitialized);
    symm(side, uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto symm(char side, char uplo, T&& A, U&& B)
{
    return symm(side, uplo, 1.0, A, B);
}

/**
 * Perform the matrix multiplication \f$ C = \alpha AB + \beta C \f$ or
 * \f$ C = \alpha BA + \beta C \f$, where A is a Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
hemm(char side, char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = C_.length(0);
    auto n = C_.length(1);

    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == n);
    MARRAY_ASSERT(A_.length(0) == side == 'L' ? m : n);
    MARRAY_ASSERT(A_.length(1) == side == 'L' ? m : n);

#ifdef BLIS_H
    MARRAY_ASSERT(side == 'L' || side == 'R');
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_hemm(side == 'L' ? BLIS_LEFT : BLIS_RIGHT, uplo == 'U' ? BLIS_UPPER : BLIS_LOWER,
             BLIS_NO_CONJUGATE, BLIS_NO_TRANSPOSE, m, n,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    B_.data(), B_.stride(0), B_.stride(1),
              beta, C_.data(), C_.stride(0), C_.stride(1));
#else
    if (C.stride(0) > 1)
    {
        side = side == 'L' ? 'R' :
               side == 'R' ? 'L' : side;
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        std::swap(m, n);
        A_.transpose();
        B_.transpose();
        C_.transpose();
    }

    MARRAY_ASSERT(A_.stride(0) == 1 && A_.stride(1) > 1);
    MARRAY_ASSERT(B_.stride(0) == 1 && B_.stride(1) > 1);
    MARRAY_ASSERT(C_.stride(0) == 1 && C_.stride(1) > 1);

    hemm(side, uplo, m, n,
         alpha, A_.data(), A_.stride(1),
                B_.data(), B_.stride(1),
          beta, C_.data(), C_.stride(1));
#endif
}

/**
 * Perform the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
hemm(char side, char uplo, T&& A, U&& B, V&& C)
{
    hemm(side, uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the matrix multiplication \f$ C = \alpha AB \f$ or
 * \f$ C = \alpha BA \f$, where A is a Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto hemm(char side, char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({B.length(0), B.length(1)}, uninitialized);
    hemm(side, uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto hemm(char side, char uplo, T&& A, U&& B)
{
    return hemm(side, uplo, 1.0, A, B);
}

/**
 * Perform the matrix multiplication \f$ B = \alpha AB \f$ or
 * \f$ B = \alpha BA \f$, where A is a triangular matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is lower-triangular or 'U' if A is upper-triangular.
 *
 * @param diag  UNIT if A is unit-diagonal or NONUNIT if A is general triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
trmm(char side, char uplo, char diag, T alpha, U&& A, V&& B)
{
    auto A_ = A.view();
    auto B_ = B.view();

    auto m = B_.length(0);
    auto n = B_.length(1);

    MARRAY_ASSERT(A_.length(0) == side == 'L' ? m : n);
    MARRAY_ASSERT(A_.length(1) == side == 'L' ? m : n);

#ifdef BLIS_H
    MARRAY_ASSERT(side == 'L' || side == 'R');
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    MARRAY_ASSERT(diag == 'N' || diag == 'U');
    bli_trmm(side == 'L' ? BLIS_LEFT : BLIS_RIGHT, uplo == 'U' ? BLIS_UPPER : BLIS_LOWER,
             BLIS_NO_TRANSPOSE, diag == 'N' ? BLIS_NONUNIT_DIAG : BLIS_UNIT_DIAG, m, n,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    B_.data(), B_.stride(0), B_.stride(1));
#else
    if (A.stride(0) > 1)
    {
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        A_.transpose();
    }

    if (B.stride(0) > 1)
    {
        side = side == 'L' ? 'R' :
               side == 'R' ? 'L' : side;
        std::swap(m, n);
        B_.transpose();
    }

    MARRAY_ASSERT(A_.stride(0) == 1 && A_.stride(1) > 1);
    MARRAY_ASSERT(B_.stride(0) == 1 && B_.stride(1) > 1);

    trmm(side, uplo, diag, m, n,
         alpha, A_.data(), A_.stride(1),
                B_.data(), B_.stride(1));
#endif
}

/**
 * Solve the system of equations \f$ \alpha B = AX \f$ or
 * \f$ \alpha B = XA \f$, where A is a triangular matrix.
 *
 * @param side  'L' if `B = AX` is to be solved, or 'R' is `B = XA' is to be solved.
 *
 * @param uplo  'L' if A is lower-triangular or 'U' if A is upper-triangular.
 *
 * @param diag  UNIT if A is unit-diagonal or NONUNIT if A is general triangular.
 *
 * @param alpha Scalar factor to apply to the solution matrix.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
trsm(char side, char uplo, char diag, T alpha, U&& A, V&& B)
{
    auto A_ = A.view();
    auto B_ = B.view();

    auto m = B_.length(0);
    auto n = B_.length(1);

    MARRAY_ASSERT(A_.length(0) == side == 'L' ? m : n);
    MARRAY_ASSERT(A_.length(1) == side == 'L' ? m : n);

#ifdef BLIS_H
    MARRAY_ASSERT(side == 'L' || side == 'R');
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    MARRAY_ASSERT(diag == 'N' || diag == 'U');
    bli_trsm(side == 'L' ? BLIS_LEFT : BLIS_RIGHT, uplo == 'U' ? BLIS_UPPER : BLIS_LOWER,
             BLIS_NO_TRANSPOSE, diag == 'N' ? BLIS_NONUNIT_DIAG : BLIS_UNIT_DIAG, m, n,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    B_.data(), B_.stride(0), B_.stride(1));
#else
    if (A.stride(0) > 1)
    {
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        A_.transpose();
    }

    if (B.stride(0) > 1)
    {
        side = side == 'L' ? 'R' :
               side == 'R' ? 'L' : side;
        std::swap(m, n);
        B_.transpose();
    }

    MARRAY_ASSERT(A_.stride(0) == 1 && A_.stride(1) > 1);
    MARRAY_ASSERT(B_.stride(0) == 1 && B_.stride(1) > 1);

    trsm(side, uplo, diag, m, n,
         alpha, A_.data(), A_.stride(1),
                B_.data(), B_.stride(1));
#endif
}

/**
 * Perform the symmetric matrix multiplication \f$ C = \alpha AA^T + \beta C \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AA^T`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<W,2>>
syrk(char uplo, T alpha, U&& A, V beta, W&& C)
{
    auto A_ = A.view();
    auto C_ = C.view();

    auto m = A_.length(0);
    auto k = A_.length(1);

    MARRAY_ASSERT(C_.length(0) == m);
    MARRAY_ASSERT(C_.length(1) == m);

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_syrk(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_TRANSPOSE, m, k,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
              beta, C_.data(), C_.stride(0), C_.stride(1));
#else
    if (C_.stride(0) > 1)
    {
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        A_.transpose();
        C_.transpose();
    }

    char trans = 'N';

    if (A_.stride(0) > 1)
    {
        trans = 'T';
        A.transpose();
    }

    MARRAY_ASSERT(A_.stride(0) == 1);
    MARRAY_ASSERT(C_.stride(0) == 1);

    syrk(uplo, trans, m, k,
         alpha, A_.data(), A_.stride(1),
          beta, C_.data(), C_.stride(1));
#endif
}

/**
 * Perform the symmetric matrix multiplication \f$ C = AA^T \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2>>
syrk(char uplo, T&& A, U&& C)
{
    syrk(uplo, 1.0, A, 0.0, C);
}

/**
 * Return the result of the symmetric matrix multiplication \f$ \alpha AA^T \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AA^T`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the scaled product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2>>>
auto syrk(char uplo, T alpha, U&& A)
{
    marray<detail::value_type<U>,2> C({A.length(0), A.length(0)}, uninitialized);
    syrk(uplo, alpha, A, 0.0, C);
    return C;
}

/**
 * Return the result of the symmetric matrix multiplication \f$ AA^T \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the product.
 */
template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2>>>
auto syrk(char uplo, T&& A)
{
    return syrk(uplo, 1.0, A);
}

/**
 * Perform the symmetric matrix multiplication \f$ C = \alpha (AB^T + BA^T) + \beta C \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^T + BA^T`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
syr2k(char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = A_.length(0);
    auto k = A_.length(1);

    MARRAY_ASSERT(C_.length(0) == m);
    MARRAY_ASSERT(C_.length(1) == m);
    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == k);

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_syr2k(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, m, k,
              alpha, A_.data(), A_.stride(0), A_.stride(1),
                     B_.data(), B_.stride(0), B_.stride(1),
               beta, C_.data(), C_.stride(0), C_.stride(1));
#else
    if (C_.stride(0) > 1)
    {
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        A_.transpose();
        B_.transpose();
        C_.transpose();
        A_.swap(B_);
    }

    char trans = 'N';

    if (A_.stride(0) > 1)
    {
        trans = 'T';
        A.transpose();
        B.transpose();
    }

    MARRAY_ASSERT(A_.stride(0) == 1);
    MARRAY_ASSERT(B_.stride(0) == 1);
    MARRAY_ASSERT(C_.stride(0) == 1);

    syr2k(uplo, trans, m, k,
          alpha, A_.data(), A_.stride(1),
                 B_.data(), B_.stride(1),
           beta, C_.data(), C_.stride(1));
#endif
}

/**
 * Perform the symmetric matrix multiplication \f$ C = AB^T + BA^T \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
syr2k(char uplo, T&& A, U&& B, V&& C)
{
    syr2k(uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the symmetric matrix multiplication \f$ \alpha (AB^T + BA^T) \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^T + BA^T`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto syr2k(char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({A.length(0), B.length(0)}, uninitialized);
    syr2k(uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the symmetric matrix multiplication \f$ AB^T + BA^T \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto syr2k(char uplo, T&& A, U&& B)
{
    return syr2k(uplo, 1.0, A, B);
}

/**
 * Perform the Hermitian matrix multiplication \f$ C = \alpha AA^H + \beta C \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AA^H`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<W,2>>
herk(char uplo, T alpha, U&& A, V beta, W&& C)
{
    auto A_ = A.view();
    auto C_ = C.view();

    auto m = A_.length(0);
    auto k = A_.length(1);

    MARRAY_ASSERT(C_.length(0) == m);
    MARRAY_ASSERT(C_.length(1) == m);

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_herk(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_TRANSPOSE, m, k,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
              beta, C_.data(), C_.stride(0), C_.stride(1));
#else
    char trans = 'N';

    if (C_.stride(0) > 1)
    {
        trans = 'T';
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        A_.transpose();
        C_.transpose();
    }

    MARRAY_ASSERT(A_.stride(0) == 1);
    MARRAY_ASSERT(C_.stride(0) == 1);

    herk(uplo, trans, m, k,
         alpha, A_.data(), A_.stride(1),
          beta, C_.data(), C_.stride(1));
#endif
}

/**
 * Perform the Hermitian matrix multiplication \f$ C = AA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2>>
herk(char uplo, T&& A, U&& C)
{
    herk(uplo, 1.0, A, 0.0, C);
}

/**
 * Return the result of the Hermitian matrix multiplication \f$ \alpha AA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AA^H`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the scaled product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2>>>
auto herk(char uplo, T alpha, U&& A)
{
    marray<detail::value_type<U>,2> C({A.length(0), A.length(0)}, uninitialized);
    herk(uplo, alpha, A, 0.0, C);
    return C;
}

/**
 * Return the result of the Hermitian matrix multiplication \f$ AA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the product.
 */
template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2>>>
auto herk(char uplo, T&& A)
{
    return herk(uplo, 1.0, A);
}

/**
 * Perform the Hermitian matrix multiplication \f$ C = \alpha AB^H + \alpha' BA^H + \beta C \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^H`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
her2k(char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = A_.length(0);
    auto k = A_.length(1);

    MARRAY_ASSERT(C_.length(0) == m);
    MARRAY_ASSERT(C_.length(1) == m);
    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == k);

#ifdef BLIS_H
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_her2k(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, m, k,
              alpha, A_.data(), A_.stride(0), A_.stride(1),
                     B_.data(), B_.stride(0), B_.stride(1),
               beta, C_.data(), C_.stride(0), C_.stride(1));
#else
    char trans = 'N';

    if (C_.stride(0) > 1)
    {
        trans = 'T';
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;
        A_.transpose();
        B_.transpose();
        C_.transpose();
        A_.swap(B_);
    }

    MARRAY_ASSERT(A_.stride(0) == 1);
    MARRAY_ASSERT(B_.stride(0) == 1);
    MARRAY_ASSERT(C_.stride(0) == 1);

    her2k(uplo, trans, m, k,
          alpha, A_.data(), A_.stride(1),
                 B_.data(), B_.stride(1),
           beta, C_.data(), C_.stride(1));
#endif
}

/**
 * Perform the Hermitian matrix multiplication \f$ C = AB^H + BA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
her2k(char uplo, T&& A, U&& B, V&& C)
{
    her2k(uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the Hermitian matrix multiplication \f$ \alpha AB^H + \alpha' BA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^H`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto her2k(char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({A.length(0), B.length(0)}, uninitialized);
    her2k(uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the Hermitian matrix multiplication \f$ AB^H + BA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto her2k(char uplo, T&& A, U&& B)
{
    return her2k(uplo, 1.0, A, B);
}

#ifdef BLIS_H

/**
 * Perform the matrix multiplication \f$ C = \alpha AB + \beta C \f$ or
 * \f$ C = \alpha BA + \beta C \f$, where A is a skew-symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
skmm(char side, char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = C_.length(0);
    auto n = C_.length(1);

    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == n);
    MARRAY_ASSERT(A_.length(0) == side == 'L' ? m : n);
    MARRAY_ASSERT(A_.length(1) == side == 'L' ? m : n);

    MARRAY_ASSERT(side == 'L' || side == 'R');
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_skmm(side == 'L' ? BLIS_LEFT : BLIS_RIGHT, uplo == 'U' ? BLIS_UPPER : BLIS_LOWER,
             BLIS_NO_CONJUGATE, BLIS_NO_TRANSPOSE, m, n,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    B_.data(), B_.stride(0), B_.stride(1),
              beta, C_.data(), C_.stride(0), C_.stride(1));
}

/**
 * Perform the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a skew-symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
skmm(char side, char uplo, T&& A, U&& B, V&& C)
{
    skmm(side, uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the matrix multiplication \f$ C = \alpha AB \f$ or
 * \f$ C = \alpha BA \f$, where A is a skew-symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto skmm(char side, char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({B.length(0), B.length(1)}, uninitialized);
    skmm(side, uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a skew-symmetric matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto skmm(char side, char uplo, T&& A, U&& B)
{
    return skmm(side, uplo, 1.0, A, B);
}

/**
 * Perform the matrix multiplication \f$ C = \alpha AB + \beta C \f$ or
 * \f$ C = \alpha BA + \beta C \f$, where A is a skew-Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
shmm(char side, char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = C_.length(0);
    auto n = C_.length(1);

    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == n);
    MARRAY_ASSERT(A_.length(0) == side == 'L' ? m : n);
    MARRAY_ASSERT(A_.length(1) == side == 'L' ? m : n);

    MARRAY_ASSERT(side == 'L' || side == 'R');
    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_shmm(side == 'L' ? BLIS_LEFT : BLIS_RIGHT, uplo == 'U' ? BLIS_UPPER : BLIS_LOWER,
             BLIS_NO_CONJUGATE, BLIS_NO_TRANSPOSE, m, n,
             alpha, A_.data(), A_.stride(0), A_.stride(1),
                    B_.data(), B_.stride(0), B_.stride(1),
              beta, C_.data(), C_.stride(0), C_.stride(1));
}

/**
 * Perform the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a skew-Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
shmm(char side, char uplo, T&& A, U&& B, V&& C)
{
    shmm(side, uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the matrix multiplication \f$ C = \alpha AB \f$ or
 * \f$ C = \alpha BA \f$, where A is a skew-Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB` or `BA`.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto shmm(char side, char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({B.length(0), B.length(1)}, uninitialized);
    shmm(side, uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the matrix multiplication \f$ C = AB \f$ or
 * \f$ C = BA \f$, where A is a skew-Hermitian matrix.
 *
 * @param side  'L' if `AB` is to be computed, or 'R' is `BA' is to be computed.
 *
 * @param uplo  'L' if A is stored lower-triangular or 'U' if A is stored upper-triangular.
 *
 * @param A     A `m`x`m` (side == 'L') or `n`x`n` (side == 'R') matrix or matrix view.
 *              Must have either a row or column stride of one.
 *
 * @param B     A `m`x`n` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto shmm(char side, char uplo, T&& A, U&& B)
{
    return shmm(side, uplo, 1.0, A, B);
}

/**
 * Perform the upper or lower triangular portion of the matrix multiplication \f$ C = \alpha AB + \beta C \f$.
 *
 * @param uplo  'L' if C is lower-triangular or 'U' if C is upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<X,2>>
gemmt(char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = A_.length(0);
    auto k = A_.length(1);

    MARRAY_ASSERT(C_.length(0) == m);
    MARRAY_ASSERT(C_.length(1) == m);
    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == k);

    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    bli_gemmt(uplo == 'U' ? BLIS_UPPER : BLIS_LOWER, BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, m, k,
              alpha, A_.data(), A_.stride(0), A_.stride(1),
                     B_.data(), B_.stride(0), B_.stride(1),
               beta, C_.data(), C_.stride(0), C_.stride(1));
}

/**
 * Perform the upper or lower triangular portion of the matrix multiplication \f$ C = AB \f$.
 *
 * @param uplo  'L' if C is lower-triangular or 'U' if C is upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
gemmt(char uplo, T&& A, U&& B, V&& C)
{
    gemmt(uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the upper or lower triangular portion of the matrix multiplication \f$ \alpha AB \f$.
 *
 * @param uplo  'L' if C is lower-triangular or 'U' if C is upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto gemmt(char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({A.length(0), B.length(0)}, uninitialized);
    gemmt(uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the upper or lower triangular portion of the matrix multiplication \f$ AB \f$.
 *
 * @param uplo  'L' if C is lower-triangular or 'U' if C is upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto gemmt(char uplo, T&& A, U&& B)
{
    return gemmt(uplo, 1.0, A, B);
}

/**
 * Perform the skew-symmetric matrix multiplication \f$ C = \alpha (AB^T - BA^T) + \beta C \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^T`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
skr2k(char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = A_.length(0);
    auto k = A_.length(1);

    MARRAY_ASSERT(C_.length(0) == m);
    MARRAY_ASSERT(C_.length(1) == m);
    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == k);

    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    //TODO
}

/**
 * Perform the skew-symmetric matrix multiplication \f$ C = AB^T - BA^T \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
skr2k(char uplo, T&& A, U&& B, V&& C)
{
    skr2k(uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the skew-symmetric matrix multiplication \f$ \alpha (AB^T - BA^T) \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^T`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto skr2k(char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({A.length(0), B.length(0)}, uninitialized);
    skr2k(uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the skew-symmetric matrix multiplication \f$ AB^T - BA^T \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto skr2k(char uplo, T&& A, U&& B)
{
    return skr2k(uplo, 1.0, A, B);
}

/**
 * Perform the skew-Hermitian matrix multiplication \f$ C = \alpha AB^H - \alpha' BA^H + \beta C \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^H`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param beta  Scalar factor for the original matrix `C`.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V, typename W, typename X>
std::enable_if_t<detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2> &&
                 detail::is_marray_like_v<W,2>>
shr2k(char uplo, T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = A_.length(0);
    auto k = A_.length(1);

    MARRAY_ASSERT(C_.length(0) == m);
    MARRAY_ASSERT(C_.length(1) == m);
    MARRAY_ASSERT(B_.length(0) == m);
    MARRAY_ASSERT(B_.length(1) == k);

    MARRAY_ASSERT(uplo == 'L' || uplo == 'U');
    //TODO
}

/**
 * Perform the skew-Hermitian matrix multiplication \f$ C = AB^H - BA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param C     A `m`x`m` matrix or matrix view. Must have either a row or
 *              column stride of one.
 */
template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<T,2> &&
                 detail::is_marray_like_v<U,2> &&
                 detail::is_marray_like_v<V,2>>
shr2k(char uplo, T&& A, U&& B, V&& C)
{
    shr2k(uplo, 1.0, A, B, 0.0, C);
}

/**
 * Return the result of the skew-Hermitian matrix multiplication \f$ \alpha AB^H - \alpha' BA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param alpha Scalar factor for the product `AB^H`.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the scaled product.
 */
template <typename T, typename U, typename V, typename=
    std::enable_if_t<detail::is_marray_like_v<U,2> &&
                     detail::is_marray_like_v<V,2>>>
auto shr2k(char uplo, T alpha, U&& A, V&& B)
{
    marray<detail::value_type<U>,2> C({A.length(0), B.length(0)}, uninitialized);
    shr2k(uplo, alpha, A, B, 0.0, C);
    return C;
}

/**
 * Return the result of the skew-Hermitian matrix multiplication \f$ AB^H - BA^H \f$.
 *
 * @param uplo  'L' if C is stored lower-triangular or 'U' if C is stored upper-triangular.
 *
 * @param A     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @param B     A `m`x`k` matrix or matrix view. Must have either a row or
 *              column stride of one.
 *
 * @return      A `m`x`m` matrix holding the product.
 */
template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,2> &&
                     detail::is_marray_like_v<U,2>>>
auto shr2k(char uplo, T&& A, U&& B)
{
    return shr2k(uplo, 1.0, A, B);
}

#endif

} //namespace blas
} //namespace MArray

#endif //__cplusplus

#endif //MARRAY_BLAS_H
