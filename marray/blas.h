#ifndef MARRAY_BLAS_HPP
#define MARRAY_BLAS_HPP

#include "detail/fortran.h"

#ifdef __cplusplus
extern "C"
{
#endif

/******************************************************************************
 *
 * Level 1 BLAS, FORTRAN prototypes
 *
 *****************************************************************************/
void    MARRAY_FC_FUNC(srotg,SROTG)  (float* a, float* b, float* c, float* s);
void    MARRAY_FC_FUNC(srotmg,SROTMG)(float* d1, float* d2, float* a, const float* b, float* param);
void    MARRAY_FC_FUNC(srot,SROT)    (const integer* n,                           float* x, const integer* incx,       float* y, const integer* incy, const float* c, const float* s);
void    MARRAY_FC_FUNC(srotm,SROTM)  (const integer* n,                           float* x, const integer* incx,       float* y, const integer* incy, float* param);
void    MARRAY_FC_FUNC(sswap,SSWAP)  (const integer* n,                           float* x, const integer* incx,       float* y, const integer* incy);
void    MARRAY_FC_FUNC(sscal,SSCAL)  (const integer* n, const float* alpha,       float* x, const integer* incx);
void    MARRAY_FC_FUNC(scopy,SCOPY)  (const integer* n,                     const float* x, const integer* incx,       float* y, const integer* incy);
void    MARRAY_FC_FUNC(saxpy,SAXPY)  (const integer* n, const float* alpha, const float* x, const integer* incx,       float* y, const integer* incy);
float   MARRAY_FC_FUNC(sdot,SDOT)    (const integer* n,                     const float* x, const integer* incx, const float* y, const integer* incy);
float   MARRAY_FC_FUNC(snrm2,SNRM2)  (const integer* n,                     const float* x, const integer* incx);
float   MARRAY_FC_FUNC(sasum,SASUM)  (const integer* n,                     const float* x, const integer* incx);
integer MARRAY_FC_FUNC(isamax,ISAMAX)(const integer* n,                     const float* x, const integer* incx);

void    MARRAY_FC_FUNC(drotg,DROTG)  (double* a, double* b, double* c, double* s);
void    MARRAY_FC_FUNC(drotmg,DROTMG)(double* d1, double* d2, double* a, const double* b, double* param);
void    MARRAY_FC_FUNC(drot,DROT)    (const integer* n,                            double* x, const integer* incx,       double* y, const integer* incy, const double* c, const double* s);
void    MARRAY_FC_FUNC(drotm,DROTM)  (const integer* n,                            double* x, const integer* incx,       double* y, const integer* incy, double* param);
void    MARRAY_FC_FUNC(dswap,DSWAP)  (const integer* n,                            double* x, const integer* incx,       double* y, const integer* incy);
void    MARRAY_FC_FUNC(dscal,DSCAL)  (const integer* n, const double* alpha,       double* x, const integer* incx);
void    MARRAY_FC_FUNC(dcopy,DCOPY)  (const integer* n,                      const double* x, const integer* incx,       double* y, const integer* incy);
void    MARRAY_FC_FUNC(daxpy,DAXPY)  (const integer* n, const double* alpha, const double* x, const integer* incx,       double* y, const integer* incy);
double  MARRAY_FC_FUNC(ddot,DDOT)    (const integer* n,                      const double* x, const integer* incx, const double* y, const integer* incy);
double  MARRAY_FC_FUNC(dnrm2,DNRM2)  (const integer* n,                      const double* x, const integer* incx);
double  MARRAY_FC_FUNC(dasum,DASUM)  (const integer* n,                      const double* x, const integer* incx);
integer MARRAY_FC_FUNC(idamax,IDAMAX)(const integer* n,                      const double* x, const integer* incx);

void     MARRAY_FC_FUNC(crotg,CROTG)  (scomplex* a, scomplex* b, float* c, scomplex* s);
void     MARRAY_FC_FUNC(csrot,CSROT)  (const integer* n,                              scomplex* x, const integer* incx,       scomplex* y, const integer* incy, const float* c, const float* s);
void     MARRAY_FC_FUNC(cswap,CSWAP)  (const integer* n,                              scomplex* x, const integer* incx,       scomplex* y, const integer* incy);
void     MARRAY_FC_FUNC(cscal,CSCAL)  (const integer* n, const scomplex* alpha,       scomplex* x, const integer* incx);
void     MARRAY_FC_FUNC(csscal,CSSCAL)(const integer* n, const    float* alpha,       scomplex* x, const integer* incx);
void     MARRAY_FC_FUNC(ccopy,CCOPY)  (const integer* n,                        const scomplex* x, const integer* incx,       scomplex* y, const integer* incy);
void     MARRAY_FC_FUNC(caxpy,CAXPY)  (const integer* n, const scomplex* alpha, const scomplex* x, const integer* incx,       scomplex* y, const integer* incy);
scomplex_f MARRAY_FC_FUNC(cdotu,CDOTU)  (const integer* n,                        const scomplex* x, const integer* incx, const scomplex* y, const integer* incy);
scomplex_f MARRAY_FC_FUNC(cdotc,CDOTC)  (const integer* n,                        const scomplex* x, const integer* incx, const scomplex* y, const integer* incy);
float    MARRAY_FC_FUNC(scnrm2,SCNRM2)(const integer* n,                        const scomplex* x, const integer* incx);
float    MARRAY_FC_FUNC(scasum,SCASUM)(const integer* n,                        const scomplex* x, const integer* incx);
integer  MARRAY_FC_FUNC(icamax,ICAMAX)(const integer* n,                        const scomplex* x, const integer* incx);

void     MARRAY_FC_FUNC(zrotg,ZROTG)  (dcomplex* a, dcomplex* b, double* c, dcomplex* s);
void     MARRAY_FC_FUNC(zdrot,ZDROT)  (const integer* n,                              dcomplex* x, const integer* incx,       dcomplex* y, const integer* incy, const double* c, const double* s);
void     MARRAY_FC_FUNC(zswap,ZSWAP)  (const integer* n,                              dcomplex* x, const integer* incx,       dcomplex* y, const integer* incy);
void     MARRAY_FC_FUNC(zscal,ZSCAL)  (const integer* n, const dcomplex* alpha,       dcomplex* x, const integer* incx);
void     MARRAY_FC_FUNC(zdscal,ZDSCAL)(const integer* n, const   double* alpha,       dcomplex* x, const integer* incx);
void     MARRAY_FC_FUNC(zcopy,ZCOPY)  (const integer* n,                        const dcomplex* x, const integer* incx,       dcomplex* y, const integer* incy);
void     MARRAY_FC_FUNC(zaxpy,ZAXPY)  (const integer* n, const dcomplex* alpha, const dcomplex* x, const integer* incx,       dcomplex* y, const integer* incy);
dcomplex_f MARRAY_FC_FUNC(zdotu,ZDOTU)  (const integer* n,                        const dcomplex* x, const integer* incx, const dcomplex* y, const integer* incy);
dcomplex_f MARRAY_FC_FUNC(zdotc,ZDOTC)  (const integer* n,                        const dcomplex* x, const integer* incx, const dcomplex* y, const integer* incy);
double   MARRAY_FC_FUNC(dznrm2,DZNRM2)(const integer* n,                        const dcomplex* x, const integer* incx);
double   MARRAY_FC_FUNC(dzasum,DZASUM)(const integer* n,                        const dcomplex* x, const integer* incx);
integer  MARRAY_FC_FUNC(izamax,IZAMAX)(const integer* n,                        const dcomplex* x, const integer* incx);

/******************************************************************************
 *
 * Level 2 BLAS, FORTRAN prototypes
 *
 *****************************************************************************/
void MARRAY_FC_FUNC(sgemv,SGEMV)(                  const char* trans,                   const integer* m, const integer* n,                                       const float* alpha, const float* a, const integer* lda,  const float* x, const integer* incx, const float* beta, float* y, const integer* incy);
void MARRAY_FC_FUNC(sgbmv,SGBMV)(                  const char* trans,                   const integer* m, const integer* n, const integer* kl, const integer* ku, const float* alpha, const float* a, const integer* lda,  const float* x, const integer* incx, const float* beta, float* y, const integer* incy);
void MARRAY_FC_FUNC(ssymv,SSYMV)(const char* uplo,                                                        const integer* n,                                       const float* alpha, const float* a, const integer* lda,  const float* x, const integer* incx, const float* beta, float* y, const integer* incy);
void MARRAY_FC_FUNC(ssbmv,SSBMV)(const char* uplo,                                                        const integer* n, const integer* k,                     const float* alpha, const float* a, const integer* lda,  const float* x, const integer* incx, const float* beta, float* y, const integer* incy);
void MARRAY_FC_FUNC(sspmv,SSPMV)(const char* uplo,                                                        const integer* n,                                       const float* alpha, const float* ap,                     const float* x, const integer* incx, const float* beta, float* y, const integer* incy);
void MARRAY_FC_FUNC(strmv,STRMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const float* a, const integer* lda,        float* x, const integer* incx);
void MARRAY_FC_FUNC(stbmv,STBMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                         const float* a, const integer* lda,        float* x, const integer* incx);
void MARRAY_FC_FUNC(stpmv,STPMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const float* ap,                           float* x, const integer* incx);
void MARRAY_FC_FUNC(strsv,STRSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const float* a, const integer* lda,        float* x, const integer* incx);
void MARRAY_FC_FUNC(stbsv,STBSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                         const float* a, const integer* lda,        float* x, const integer* incx);
void MARRAY_FC_FUNC(stpsv,STPSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                           const float* ap,                           float* x, const integer* incx);
void MARRAY_FC_FUNC(sger,SGER)  (                                                       const integer* m, const integer* n,                                       const float* alpha, const float* x, const integer* incx, const float* y, const integer* incy, float* a, const integer* lda);
void MARRAY_FC_FUNC(ssyr,SSYR)  (const char* uplo,                                                        const integer* n,                                       const float* alpha, const float* x, const integer* incx,                                      float* a, const integer* lda);
void MARRAY_FC_FUNC(sspr,SSPR)  (const char* uplo,                                                        const integer* n,                                       const float* alpha, const float* x, const integer* incx,                                      float* ap);
void MARRAY_FC_FUNC(ssyr2,SSYR2)(const char* uplo,                                                        const integer* n,                                       const float* alpha, const float* x, const integer* incx, const float* y, const integer* incy, float* a, const integer* lda);
void MARRAY_FC_FUNC(sspr2,SSPR2)(const char* uplo,                                                        const integer* n,                                       const float* alpha, const float* x, const integer* incx, const float* y, const integer* incy, float* ap);

void MARRAY_FC_FUNC(dgemv,DGEMV)(                  const char* trans,                   const integer* m, const integer* n,                                       const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void MARRAY_FC_FUNC(dgbmv,DGBMV)(                  const char* trans,                   const integer* m, const integer* n, const integer* kl, const integer* ku, const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void MARRAY_FC_FUNC(dsymv,DSYMV)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void MARRAY_FC_FUNC(dsbmv,DSBMV)(const char* uplo,                                                        const integer* n, const integer* k,                     const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void MARRAY_FC_FUNC(dspmv,DSPMV)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* ap,                    const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void MARRAY_FC_FUNC(dtrmv,DTRMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* a, const integer* lda,       double* x, const integer* incx);
void MARRAY_FC_FUNC(dtbmv,DTBMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                          const double* a, const integer* lda,       double* x, const integer* incx);
void MARRAY_FC_FUNC(dtpmv,DTPMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* ap,                          double* x, const integer* incx);
void MARRAY_FC_FUNC(dtrsv,DTRSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* a, const integer* lda,       double* x, const integer* incx);
void MARRAY_FC_FUNC(dtbsv,DTBSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                          const double* a, const integer* lda,       double* x, const integer* incx);
void MARRAY_FC_FUNC(dtpsv,DTPSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* ap,                          double* x, const integer* incx);
void MARRAY_FC_FUNC(dger,DGER)  (                                                       const integer* m, const integer* n,                                       const double* alpha, const double* x, const integer* incx, const double* y, const integer* incy, double* a, const integer* lda);
void MARRAY_FC_FUNC(dsyr,DSYR)  (const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx,                                       double* a, const integer* lda);
void MARRAY_FC_FUNC(dspr,DSPR)  (const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx,                                       double* ap);
void MARRAY_FC_FUNC(dsyr2,DSYR2)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx, const double* y, const integer* incy, double* a, const integer* lda);
void MARRAY_FC_FUNC(dspr2,DSPR2)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx, const double* y, const integer* incy, double* ap);

void MARRAY_FC_FUNC(cgemv,CGEMV)(                  const char* trans,                   const integer* m, const integer* n,                                       const scomplex* alpha, const scomplex* a, const integer* lda,  const scomplex* x, const integer* incx, const scomplex* beta, scomplex* y, const integer* incy);
void MARRAY_FC_FUNC(cgbmv,CGBMV)(                  const char* trans,                   const integer* m, const integer* n, const integer* kl, const integer* ku, const scomplex* alpha, const scomplex* a, const integer* lda,  const scomplex* x, const integer* incx, const scomplex* beta, scomplex* y, const integer* incy);
void MARRAY_FC_FUNC(chemv,CHEMV)(const char* uplo,                                                        const integer* n,                                       const scomplex* alpha, const scomplex* a, const integer* lda,  const scomplex* x, const integer* incx, const scomplex* beta, scomplex* y, const integer* incy);
void MARRAY_FC_FUNC(chbmv,CHBMV)(const char* uplo,                                                        const integer* n, const integer* k,                     const scomplex* alpha, const scomplex* a, const integer* lda,  const scomplex* x, const integer* incx, const scomplex* beta, scomplex* y, const integer* incy);
void MARRAY_FC_FUNC(chpmv,CHPMV)(const char* uplo,                                                        const integer* n,                                       const scomplex* alpha, const scomplex* ap,                     const scomplex* x, const integer* incx, const scomplex* beta, scomplex* y, const integer* incy);
void MARRAY_FC_FUNC(ctrmv,CTRMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const scomplex* a, const integer* lda,        scomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ctbmv,CTBMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                            const scomplex* a, const integer* lda,        scomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ctpmv,CTPMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const scomplex* ap,                           scomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ctrsv,CTRSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const scomplex* a, const integer* lda,        scomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ctbsv,CTBSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                            const scomplex* a, const integer* lda,        scomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ctpsv,CTPSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const scomplex* ap,                           scomplex* x, const integer* incx);
void MARRAY_FC_FUNC(cgeru,CGERU)(                                                       const integer* m, const integer* n,                                       const scomplex* alpha, const scomplex* x, const integer* incx, const scomplex* y, const integer* incy, scomplex* a, const integer* lda);
void MARRAY_FC_FUNC(cgerc,CGERC)(                                                       const integer* m, const integer* n,                                       const scomplex* alpha, const scomplex* x, const integer* incx, const scomplex* y, const integer* incy, scomplex* a, const integer* lda);
void MARRAY_FC_FUNC(cher,CHER)  (const char* uplo,                                                        const integer* n,                                       const    float* alpha, const scomplex* x, const integer* incx,                                         scomplex* a, const integer* lda);
void MARRAY_FC_FUNC(chpr,CHPR)  (const char* uplo,                                                        const integer* n,                                       const    float* alpha, const scomplex* x, const integer* incx,                                         scomplex* ap);
void MARRAY_FC_FUNC(cher2,CHER2)(const char* uplo,                                                        const integer* n,                                       const scomplex* alpha, const scomplex* x, const integer* incx, const scomplex* y, const integer* incy, scomplex* a, const integer* lda);
void MARRAY_FC_FUNC(chpr2,CHPR2)(const char* uplo,                                                        const integer* n,                                       const scomplex* alpha, const scomplex* x, const integer* incx, const scomplex* y, const integer* incy, scomplex* ap);

void MARRAY_FC_FUNC(zgemv,ZGEMV)(                  const char* trans,                   const integer* m, const integer* n,                                       const dcomplex* alpha, const dcomplex* a, const integer* lda,  const dcomplex* x, const integer* incx, const dcomplex* beta, dcomplex* y, const integer* incy);
void MARRAY_FC_FUNC(zgbmv,ZGBMV)(                  const char* trans,                   const integer* m, const integer* n, const integer* kl, const integer* ku, const dcomplex* alpha, const dcomplex* a, const integer* lda,  const dcomplex* x, const integer* incx, const dcomplex* beta, dcomplex* y, const integer* incy);
void MARRAY_FC_FUNC(zhemv,ZHEMV)(const char* uplo,                                                        const integer* n,                                       const dcomplex* alpha, const dcomplex* a, const integer* lda,  const dcomplex* x, const integer* incx, const dcomplex* beta, dcomplex* y, const integer* incy);
void MARRAY_FC_FUNC(zhbmv,ZHBMV)(const char* uplo,                                                        const integer* n, const integer* k,                     const dcomplex* alpha, const dcomplex* a, const integer* lda,  const dcomplex* x, const integer* incx, const dcomplex* beta, dcomplex* y, const integer* incy);
void MARRAY_FC_FUNC(zhpmv,ZHPMV)(const char* uplo,                                                        const integer* n,                                       const dcomplex* alpha, const dcomplex* ap,                     const dcomplex* x, const integer* incx, const dcomplex* beta, dcomplex* y, const integer* incy);
void MARRAY_FC_FUNC(ztrmv,ZTRMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const dcomplex* a, const integer* lda,        dcomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ztbmv,ZTBMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                            const dcomplex* a, const integer* lda,        dcomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ztpmv,ZTPMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const dcomplex* ap,                           dcomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ztrsv,ZTRSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const dcomplex* a, const integer* lda,        dcomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ztbsv,ZTBSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                            const dcomplex* a, const integer* lda,        dcomplex* x, const integer* incx);
void MARRAY_FC_FUNC(ztpsv,ZTPSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                              const dcomplex* ap,                           dcomplex* x, const integer* incx);
void MARRAY_FC_FUNC(zgerc,ZGERC)(                                                       const integer* m, const integer* n,                                       const dcomplex* alpha, const dcomplex* x, const integer* incx, const dcomplex* y, const integer* incy, dcomplex* a, const integer* lda);
void MARRAY_FC_FUNC(zgeru,ZGERU)(                                                       const integer* m, const integer* n,                                       const dcomplex* alpha, const dcomplex* x, const integer* incx, const dcomplex* y, const integer* incy, dcomplex* a, const integer* lda);
void MARRAY_FC_FUNC(zher,ZHER)  (const char* uplo,                                                        const integer* n,                                       const   double* alpha, const dcomplex* x, const integer* incx,                                         dcomplex* a, const integer* lda);
void MARRAY_FC_FUNC(zhpr,ZHPR)  (const char* uplo,                                                        const integer* n,                                       const   double* alpha, const dcomplex* x, const integer* incx,                                         dcomplex* ap);
void MARRAY_FC_FUNC(zher2,ZHER2)(const char* uplo,                                                        const integer* n,                                       const dcomplex* alpha, const dcomplex* x, const integer* incx, const dcomplex* y, const integer* incy, dcomplex* a, const integer* lda);
void MARRAY_FC_FUNC(zhpr2,ZHPR2)(const char* uplo,                                                        const integer* n,                                       const dcomplex* alpha, const dcomplex* x, const integer* incx, const dcomplex* y, const integer* incy, dcomplex* ap);

/******************************************************************************
 *
 * Level 3 BLAS, FORTRAN prototypes
 *
 *****************************************************************************/
void MARRAY_FC_FUNC(sgemm,SGEMM)  (                                    const char* transa, const char* transb,                   const integer* m, const integer* n, const integer* k, const float* alpha, const float* a, const integer* lda, const float* b, const integer* ldb, const float* beta, float* c, const integer* ldc);
void MARRAY_FC_FUNC(ssymm,SSYMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const float* alpha, const float* a, const integer* lda, const float* b, const integer* ldb, const float* beta, float* c, const integer* ldc);
void MARRAY_FC_FUNC(ssyrk,SSYRK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const float* alpha, const float* a, const integer* lda,                                     const float* beta, float* c, const integer* ldc);
void MARRAY_FC_FUNC(ssyr2k,SSYR2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const float* alpha, const float* a, const integer* lda, const float* b, const integer* ldb, const float* beta, float* c, const integer* ldc);
void MARRAY_FC_FUNC(strmm,STRMM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const float* alpha, const float* a, const integer* lda,       float* b, const integer* ldb);
void MARRAY_FC_FUNC(strsm,STRSM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const float* alpha, const float* a, const integer* lda,       float* b, const integer* ldb);

void MARRAY_FC_FUNC(dgemm,DGEMM)  (                                    const char* transa, const char* transb,                   const integer* m, const integer* n, const integer* k, const double* alpha, const double* a, const integer* lda, const double* b, const integer* ldb, const double* beta, double* c, const integer* ldc);
void MARRAY_FC_FUNC(dsymm,DSYMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const double* alpha, const double* a, const integer* lda, const double* b, const integer* ldb, const double* beta, double* c, const integer* ldc);
void MARRAY_FC_FUNC(dsyrk,DSYRK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const double* alpha, const double* a, const integer* lda,                                      const double* beta, double* c, const integer* ldc);
void MARRAY_FC_FUNC(dsyr2k,DSYR2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const double* alpha, const double* a, const integer* lda, const double* b, const integer* ldb, const double* beta, double* c, const integer* ldc);
void MARRAY_FC_FUNC(dtrmm,DTRMM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const double* alpha, const double* a, const integer* lda,       double* b, const integer* ldb);
void MARRAY_FC_FUNC(dtrsm,DTRSM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const double* alpha, const double* a, const integer* lda,       double* b, const integer* ldb);

void MARRAY_FC_FUNC(cgemm,CGEMM)  (                                    const char* transa, const char* transb,                   const integer* m, const integer* n, const integer* k, const scomplex* alpha, const scomplex* a, const integer* lda, const scomplex* b, const integer* ldb, const scomplex* beta, scomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(chemm,CHEMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const scomplex* alpha, const scomplex* a, const integer* lda, const scomplex* b, const integer* ldb, const scomplex* beta, scomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(csymm,CSYMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const scomplex* alpha, const scomplex* a, const integer* lda, const scomplex* b, const integer* ldb, const scomplex* beta, scomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(csyrk,CSYRK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const scomplex* alpha, const scomplex* a, const integer* lda,                                        const scomplex* beta, scomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(csyr2k,CSYR2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const scomplex* alpha, const scomplex* a, const integer* lda, const scomplex* b, const integer* ldb, const scomplex* beta, scomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(cherk,CHERK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const    float* alpha, const scomplex* a, const integer* lda,                                        const    float* beta, scomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(cher2k,CHER2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const scomplex* alpha, const scomplex* a, const integer* lda, const scomplex* b, const integer* ldb, const    float* beta, scomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(ctrmm,CTRMM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const scomplex* alpha, const scomplex* a, const integer* lda,       scomplex* b, const integer* ldb);
void MARRAY_FC_FUNC(ctrsm,CTRSM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const scomplex* alpha, const scomplex* a, const integer* lda,       scomplex* b, const integer* ldb);

void MARRAY_FC_FUNC(zgemm,ZGEMM)  (                                    const char* transa, const char* transb,                   const integer* m, const integer* n, const integer* k, const dcomplex* alpha, const dcomplex* a, const integer* lda, const dcomplex* b, const integer* ldb, const dcomplex* beta, dcomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(zhemm,ZHEMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const dcomplex* alpha, const dcomplex* a, const integer* lda, const dcomplex* b, const integer* ldb, const dcomplex* beta, dcomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(zsymm,ZSYMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const dcomplex* alpha, const dcomplex* a, const integer* lda, const dcomplex* b, const integer* ldb, const dcomplex* beta, dcomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(zsyrk,ZSYRK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const dcomplex* alpha, const dcomplex* a, const integer* lda,                                        const dcomplex* beta, dcomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(zsyr2k,ZSYR2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const dcomplex* alpha, const dcomplex* a, const integer* lda, const dcomplex* b, const integer* ldb, const dcomplex* beta, dcomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(zherk,ZHERK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const   double* alpha, const dcomplex* a, const integer* lda,                                        const   double* beta, dcomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(zher2k,ZHER2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const dcomplex* alpha, const dcomplex* a, const integer* lda, const dcomplex* b, const integer* ldb, const   double* beta, dcomplex* c, const integer* ldc);
void MARRAY_FC_FUNC(ztrmm,ZTRMM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const dcomplex* alpha, const dcomplex* a, const integer* lda,       dcomplex* b, const integer* ldb);
void MARRAY_FC_FUNC(ztrsm,ZTRSM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const dcomplex* alpha, const dcomplex* a, const integer* lda,       dcomplex* b, const integer* ldb);

/******************************************************************************
 *
 * Level 1 BLAS, C wrappers
 *
 *****************************************************************************/

static inline void   c_srotg (float* a, float* b, float* c, float* s)
{
    MARRAY_FC_FUNC(srotg,SROTG)(a, b, c, s);
}

static inline void   c_srotmg(float* d1, float* d2, float* a, const float b, float* param)
{
    MARRAY_FC_FUNC(srotmg,SROTMG)(d1, d2, a, &b, param);
}

static inline void   c_srot  (const integer n, float* x, const integer incx,
                                               float* y, const integer incy, const float c, const float s)
{
    MARRAY_FC_FUNC(srot,SROT)(&n, x, &incx, y, &incy, &c, &s);
}

static inline void   c_srotm (const integer n, float* x, const integer incx,
                                               float* y, const integer incy, float* param)
{
    MARRAY_FC_FUNC(srotm,SROTM)(&n, x, &incx, y, &incy, param);
}

static inline void   c_sswap (const integer n, float* x, const integer incx,
                                               float* y, const integer incy)
{
    MARRAY_FC_FUNC(sswap,SSWAP)(&n, x, &incx, y, &incy);
}

static inline void   c_sscal (const integer n, const float alpha, float* x, const integer incx)
{
    MARRAY_FC_FUNC(sscal,SSCAL)(&n, &alpha, x, &incx);
}

static inline void   c_scopy (const integer n, const float* x, const integer incx,
                                                     float* y, const integer incy)
{
    MARRAY_FC_FUNC(scopy,SCOPY)(&n, x, &incx, y, &incy);
}

static inline void   c_saxpy (const integer n, const float alpha, const float* x, const integer incx,
                                                                        float* y, const integer incy)
{
    MARRAY_FC_FUNC(saxpy,SAXPY)(&n, &alpha, x, &incx, y, &incy);
}

static inline float c_sdot  (const integer n, const float* x, const integer incx,
                                              const float* y, const integer incy)
{
    return MARRAY_FC_FUNC(sdot,SDOT)(&n, x, &incx, y, &incy);
}

static inline float c_snrm2 (const integer n, const float* x, const integer incx)
{
    return MARRAY_FC_FUNC(snrm2,SNRM2)(&n, x, &incx);
}

static inline float c_sasum (const integer n, const float* x, const integer incx)
{
    return MARRAY_FC_FUNC(sasum,SASUM)(&n, x, &incx);
}

static inline integer c_isamax(const integer n, const float* x, const integer incx)
{
    return MARRAY_FC_FUNC(isamax,ISAMAX)(&n, x, &incx)-1;
}

static inline void   c_drotg (double* a, double* b, double* c, double* s)
{
    MARRAY_FC_FUNC(drotg,DROTG)(a, b, c, s);
}

static inline void   c_drotmg(double* d1, double* d2, double* a, const double b, double* param)
{
    MARRAY_FC_FUNC(drotmg,DROTMG)(d1, d2, a, &b, param);
}

static inline void   c_drot  (const integer n, double* x, const integer incx,
                                               double* y, const integer incy, const double c, const double s)
{
    MARRAY_FC_FUNC(drot,DROT)(&n, x, &incx, y, &incy, &c, &s);
}

static inline void   c_drotm (const integer n, double* x, const integer incx,
                                               double* y, const integer incy, double* param)
{
    MARRAY_FC_FUNC(drotm,DROTM)(&n, x, &incx, y, &incy, param);
}

static inline void   c_dswap (const integer n, double* x, const integer incx,
                                               double* y, const integer incy)
{
    MARRAY_FC_FUNC(dswap,DSWAP)(&n, x, &incx, y, &incy);
}

static inline void   c_dscal (const integer n, const double alpha, double* x, const integer incx)
{
    MARRAY_FC_FUNC(dscal,DSCAL)(&n, &alpha, x, &incx);
}

static inline void   c_dcopy (const integer n, const double* x, const integer incx,
                                                     double* y, const integer incy)
{
    MARRAY_FC_FUNC(dcopy,DCOPY)(&n, x, &incx, y, &incy);
}

static inline void   c_daxpy (const integer n, const double alpha, const double* x, const integer incx,
                                                                         double* y, const integer incy)
{
    MARRAY_FC_FUNC(daxpy,DAXPY)(&n, &alpha, x, &incx, y, &incy);
}

static inline double c_ddot  (const integer n, const double* x, const integer incx,
                                               const double* y, const integer incy)
{
    return MARRAY_FC_FUNC(ddot,DDOT)(&n, x, &incx, y, &incy);
}

static inline double c_dnrm2 (const integer n, const double* x, const integer incx)
{
    return MARRAY_FC_FUNC(dnrm2,DNRM2)(&n, x, &incx);
}

static inline double c_dasum (const integer n, const double* x, const integer incx)
{
    return MARRAY_FC_FUNC(dasum,DASUM)(&n, x, &incx);
}

static inline integer c_idamax(const integer n, const double* x, const integer incx)
{
    return MARRAY_FC_FUNC(idamax,IDAMAX)(&n, x, &incx)-1;
}

static inline void   c_crotg (scomplex* a, scomplex* b, float* c, scomplex* s)
{
    MARRAY_FC_FUNC(crotg,CROTG)(a, b, c, s);
}

static inline void   c_csrot  (const integer n, scomplex* x, const integer incx,
                                                scomplex* y, const integer incy, const float c, const float s)
{
    MARRAY_FC_FUNC(csrot,CSROT)(&n, x, &incx, y, &incy, &c, &s);
}

static inline void   c_cswap (const integer n, scomplex* x, const integer incx,
                                               scomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(cswap,CSWAP)(&n, x, &incx, y, &incy);
}

static inline void   c_cscal (const integer n, const scomplex alpha, scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(cscal,CSCAL)(&n, &alpha, x, &incx);
}

static inline void   c_csscal (const integer n, const float alpha, scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(csscal,CSSCAL)(&n, &alpha, x, &incx);
}

static inline void   c_ccopy (const integer n, const scomplex* x, const integer incx,
                                                     scomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(ccopy,CCOPY)(&n, x, &incx, y, &incy);
}

static inline void   c_caxpy (const integer n, const scomplex alpha, const scomplex* x, const integer incx,
                                                                           scomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(caxpy,CAXPY)(&n, &alpha, x, &incx, y, &incy);
}

static inline scomplex c_cdotu (const integer n, const scomplex* x, const integer incx,
                                                 const scomplex* y, const integer incy)
{
    scomplex_f tmp = MARRAY_FC_FUNC(cdotu,CDOTU)(&n, x, &incx, y, &incy);
    return MAKE_COMPLEX(float, tmp.real, tmp.imag);
}

static inline scomplex c_cdotc (const integer n, const scomplex* x, const integer incx,
                                                 const scomplex* y, const integer incy)
{
    scomplex_f tmp = MARRAY_FC_FUNC(cdotc,CDOTC)(&n, x, &incx, y, &incy);
    return MAKE_COMPLEX(float, tmp.real, tmp.imag);
}

static inline float c_scnrm2 (const integer n, const scomplex* x, const integer incx)
{
    return MARRAY_FC_FUNC(scnrm2,SCNRM2)(&n, x, &incx);
}

static inline float c_scasum (const integer n, const scomplex* x, const integer incx)
{
    return MARRAY_FC_FUNC(scasum,SCASUM)(&n, x, &incx);
}

static inline integer c_icamax(const integer n, const scomplex* x, const integer incx)
{
    return MARRAY_FC_FUNC(icamax,ICAMAX)(&n, x, &incx)-1;
}

static inline void   c_zrotg (dcomplex* a, dcomplex* b, double* c, dcomplex* s)
{
    MARRAY_FC_FUNC(zrotg,ZROTG)(a, b, c, s);
}

static inline void   c_zdrot  (const integer n, dcomplex* x, const integer incx,
                                                dcomplex* y, const integer incy, const double c, const double s)
{
    MARRAY_FC_FUNC(zdrot,ZDROT)(&n, x, &incx, y, &incy, &c, &s);
}

static inline void   c_zswap (const integer n, dcomplex* x, const integer incx,
                                               dcomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(zswap,ZSWAP)(&n, x, &incx, y, &incy);
}

static inline void   c_zscal (const integer n, const dcomplex alpha, dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(zscal,ZSCAL)(&n, &alpha, x, &incx);
}

static inline void   c_zdscal (const integer n, const double alpha, dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(zdscal,ZDSCAL)(&n, &alpha, x, &incx);
}

static inline void   c_zcopy (const integer n, const dcomplex* x, const integer incx,
                                                     dcomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(zcopy,ZCOPY)(&n, x, &incx, y, &incy);
}

static inline void   c_zaxpy (const integer n, const dcomplex alpha, const dcomplex* x, const integer incx,
                                                                           dcomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(zaxpy,ZAXPY)(&n, &alpha, x, &incx, y, &incy);
}

static inline dcomplex c_zdotu (const integer n, const dcomplex* x, const integer incx,
                                                 const dcomplex* y, const integer incy)
{
    dcomplex_f tmp = MARRAY_FC_FUNC(zdotu,ZDOTU)(&n, x, &incx, y, &incy);
    return MAKE_COMPLEX(double, tmp.real, tmp.imag);
}

static inline dcomplex c_zdotc (const integer n, const dcomplex* x, const integer incx,
                                                 const dcomplex* y, const integer incy)
{
    dcomplex_f tmp = MARRAY_FC_FUNC(zdotc,ZDOTC)(&n, x, &incx, y, &incy);
    return MAKE_COMPLEX(double, tmp.real, tmp.imag);
}

static inline double c_dznrm2 (const integer n, const dcomplex* x, const integer incx)
{
    return MARRAY_FC_FUNC(dznrm2,DZNRM2)(&n, x, &incx);
}

static inline double c_dzasum (const integer n, const dcomplex* x, const integer incx)
{
    return MARRAY_FC_FUNC(dzasum,DZASUM)(&n, x, &incx);
}

static inline integer c_izamax(const integer n, const dcomplex* x, const integer incx)
{
    return MARRAY_FC_FUNC(izamax,IZAMAX)(&n, x, &incx)-1;
}

/******************************************************************************
 *
 * Level 2 BLAS, C wrappers
 *
 *****************************************************************************/
static inline void c_sgemv(const char trans, const integer m, const integer n,
                           const float alpha, const float* a, const integer lda,
                                              const float* x, const integer incx,
                           const float  beta,       float* y, const integer incy)
{
    MARRAY_FC_FUNC(sgemv,SGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_sgbmv(const char trans,
                           const integer m, const integer n, const integer kl, const integer ku,
                           const float alpha, const float* a, const integer lda,
                                              const float* x, const integer incx,
                           const float  beta,       float* y, const integer incy)
{
    MARRAY_FC_FUNC(sgbmv,SGBMV)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_ssymv(const char uplo, const integer n,
                           const float alpha, const float* a, const integer lda,
                                              const float* x, const integer incx,
                           const float  beta,       float* y, const integer incy)
{
    MARRAY_FC_FUNC(ssymv,SSYMV)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_ssbmv(const char uplo, const integer n, const integer k,
                           const float alpha, const float* a, const integer lda,
                                              const float* x, const integer incx,
                           const float  beta,       float* y, const integer incy)
{
    MARRAY_FC_FUNC(ssbmv,SSBMV)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_sspmv(const char uplo, const integer n,
                           const float alpha, const float* ap,
                                              const float*  x, const integer incx,
                           const float  beta,       float*  y, const integer incy)
{
    MARRAY_FC_FUNC(sspmv,SSPMV)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
}

static inline void c_strmv(const char uplo, const char trans, const char diag, const integer n,
                           const float* a, const integer lda,
                                 float* x, const integer incx)
{
    MARRAY_FC_FUNC(strmv,STRMV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_stbmv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const float* a, const integer lda,
                                 float* x, const integer incx)
{
    MARRAY_FC_FUNC(stbmv,STBMV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_stpmv(const char uplo, const char trans, const char diag, const integer n,
                           const float* ap, float* x, const integer incx)
{
    MARRAY_FC_FUNC(stpmv,STPMV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_strsv(const char uplo, const char trans, const char diag, const integer n,
                           const float* a, const integer lda,
                                 float* x, const integer incx)
{
    MARRAY_FC_FUNC(strsv,STRSV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_stbsv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const float* a, const integer lda,
                                 float* x, const integer incx)
{
    MARRAY_FC_FUNC(stbsv,STBSV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_stpsv(const char uplo, const char trans, const char diag, const integer n,
                           const float* ap, float* x, const integer incx)
{
    MARRAY_FC_FUNC(stpsv,STPSV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_sger (const integer m, const integer n,
                           const float alpha, const float* x, const integer incx,
                                              const float* y, const integer incy,
                                                    float* a, const integer lda)
{
    MARRAY_FC_FUNC(sger,SGER)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_ssyr (const char uplo, const integer n,
                           const float alpha, const float* x, const integer incx,
                                                    float* a, const integer lda)
{
    MARRAY_FC_FUNC(ssyr,SSYR)(&uplo, &n, &alpha, x, &incx, a, &lda);
}

static inline void c_sspr (const char uplo, const integer n,
                           const float alpha, const float* x, const integer incx,
                                                    float* ap)
{
    MARRAY_FC_FUNC(sspr,SSPR)(&uplo, &n, &alpha, x, &incx, ap);
}

static inline void c_ssyr2(const char uplo, const integer n,
                           const float alpha, const float* x, const integer incx,
                                              const float* y, const integer incy,
                                                    float* a, const integer lda)
{
    MARRAY_FC_FUNC(ssyr2,SSYR2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_sspr2(const char uplo, const integer n,
                           const float alpha, const float* x, const integer incx,
                                              const float* y, const integer incy,
                                                    float* ap)
{
    MARRAY_FC_FUNC(sspr2,SSPR2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
}

static inline void c_dgemv(const char trans, const integer m, const integer n,
                           const double alpha, const double* a, const integer lda,
                                               const double* x, const integer incx,
                           const double  beta,       double* y, const integer incy)
{
    MARRAY_FC_FUNC(dgemv,DGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dgbmv(const char trans,
                           const integer m, const integer n, const integer kl, const integer ku,
                           const double alpha, const double* a, const integer lda,
                                               const double* x, const integer incx,
                           const double  beta,       double* y, const integer incy)
{
    MARRAY_FC_FUNC(dgbmv,DGBMV)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dsymv(const char uplo, const integer n,
                           const double alpha, const double* a, const integer lda,
                                               const double* x, const integer incx,
                           const double  beta,       double* y, const integer incy)
{
    MARRAY_FC_FUNC(dsymv,DSYMV)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dsbmv(const char uplo, const integer n, const integer k,
                           const double alpha, const double* a, const integer lda,
                                               const double* x, const integer incx,
                           const double  beta,       double* y, const integer incy)
{
    MARRAY_FC_FUNC(dsbmv,DSBMV)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dspmv(const char uplo, const integer n,
                           const double alpha, const double* ap,
                                               const double*  x, const integer incx,
                           const double  beta,       double*  y, const integer incy)
{
    MARRAY_FC_FUNC(dspmv,DSPMV)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
}

static inline void c_dtrmv(const char uplo, const char trans, const char diag, const integer n,
                           const double* a, const integer lda,
                                 double* x, const integer incx)
{
    MARRAY_FC_FUNC(dtrmv,DTRMV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_dtbmv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const double* a, const integer lda,
                                 double* x, const integer incx)
{
    MARRAY_FC_FUNC(dtbmv,DTBMV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_dtpmv(const char uplo, const char trans, const char diag, const integer n,
                           const double* ap, double* x, const integer incx)
{
    MARRAY_FC_FUNC(dtpmv,DTPMV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_dtrsv(const char uplo, const char trans, const char diag, const integer n,
                           const double* a, const integer lda,
                                 double* x, const integer incx)
{
    MARRAY_FC_FUNC(dtrsv,DTRSV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_dtbsv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const double* a, const integer lda,
                                 double* x, const integer incx)
{
    MARRAY_FC_FUNC(dtbsv,DTBSV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_dtpsv(const char uplo, const char trans, const char diag, const integer n,
                           const double* ap, double* x, const integer incx)
{
    MARRAY_FC_FUNC(dtpsv,DTPSV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_dger (const integer m, const integer n,
                           const double alpha, const double* x, const integer incx,
                                               const double* y, const integer incy,
                                                     double* a, const integer lda)
{
    MARRAY_FC_FUNC(dger,DGER)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_dsyr (const char uplo, const integer n,
                           const double alpha, const double* x, const integer incx,
                                                     double* a, const integer lda)
{
    MARRAY_FC_FUNC(dsyr,DSYR)(&uplo, &n, &alpha, x, &incx, a, &lda);
}

static inline void c_dspr (const char uplo, const integer n,
                           const double alpha, const double* x, const integer incx,
                                                     double* ap)
{
    MARRAY_FC_FUNC(dspr,DSPR)(&uplo, &n, &alpha, x, &incx, ap);
}

static inline void c_dsyr2(const char uplo, const integer n,
                           const double alpha, const double* x, const integer incx,
                                               const double* y, const integer incy,
                                                     double* a, const integer lda)
{
    MARRAY_FC_FUNC(dsyr2,DSYR2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_dspr2(const char uplo, const integer n,
                           const double alpha, const double* x, const integer incx,
                                               const double* y, const integer incy,
                                                     double* ap)
{
    MARRAY_FC_FUNC(dspr2,DSPR2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
}

static inline void c_cgemv(const char trans, const integer m, const integer n,
                           const scomplex alpha, const scomplex* a, const integer lda,
                                                 const scomplex* x, const integer incx,
                           const scomplex  beta,       scomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(cgemv,CGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_cgbmv(const char trans,
                           const integer m, const integer n, const integer kl, const integer ku,
                           const scomplex alpha, const scomplex* a, const integer lda,
                                                 const scomplex* x, const integer incx,
                           const scomplex  beta,       scomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(cgbmv,CGBMV)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_chemv(const char uplo, const integer n,
                           const scomplex alpha, const scomplex* a, const integer lda,
                                                 const scomplex* x, const integer incx,
                           const scomplex  beta,       scomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(chemv,CHEMV)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_chbmv(const char uplo, const integer n, const integer k,
                           const scomplex alpha, const scomplex* a, const integer lda,
                                                 const scomplex* x, const integer incx,
                           const scomplex  beta,       scomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(chbmv,CHBMV)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_chpmv(const char uplo, const integer n,
                           const scomplex alpha, const scomplex* ap,
                                                 const scomplex*  x, const integer incx,
                           const scomplex  beta,       scomplex*  y, const integer incy)
{
    MARRAY_FC_FUNC(chpmv,CHPMV)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
}

static inline void c_ctrmv(const char uplo, const char trans, const char diag, const integer n,
                           const scomplex* a, const integer lda,
                                 scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ctrmv,CTRMV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_ctbmv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const scomplex* a, const integer lda,
                                 scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ctbmv,CTBMV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_ctpmv(const char uplo, const char trans, const char diag, const integer n,
                           const scomplex* ap, scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ctpmv,CTPMV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_ctrsv(const char uplo, const char trans, const char diag, const integer n,
                           const scomplex* a, const integer lda,
                                 scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ctrsv,CTRSV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_ctbsv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const scomplex* a, const integer lda,
                                 scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ctbsv,CTBSV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_ctpsv(const char uplo, const char trans, const char diag, const integer n,
                           const scomplex* ap, scomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ctpsv,CTPSV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_cgeru(const integer m, const integer n,
                           const scomplex alpha, const scomplex* x, const integer incx,
                                                 const scomplex* y, const integer incy,
                                                       scomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(cgeru,CGERU)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_cgerc(const integer m, const integer n,
                           const scomplex alpha, const scomplex* x, const integer incx,
                                                 const scomplex* y, const integer incy,
                                                       scomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(cgerc,CGERC)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_cher (const char uplo, const integer n,
                           const float alpha, const scomplex* x, const integer incx,
                                                       scomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(cher,CHER)(&uplo, &n, &alpha, x, &incx, a, &lda);
}

static inline void c_chpr (const char uplo, const integer n,
                           const float alpha, const scomplex* x, const integer incx,
                                                       scomplex* ap)
{
    MARRAY_FC_FUNC(chpr,CHPR)(&uplo, &n, &alpha, x, &incx, ap);
}

static inline void c_cher2(const char uplo, const integer n,
                           const scomplex alpha, const scomplex* x, const integer incx,
                                                 const scomplex* y, const integer incy,
                                                       scomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(cher2,CHER2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_chpr2(const char uplo, const integer n,
                           const scomplex alpha, const scomplex* x, const integer incx,
                                                 const scomplex* y, const integer incy,
                                                       scomplex* ap)
{
    MARRAY_FC_FUNC(chpr2,CHPR2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
}

static inline void c_zgemv(const char trans, const integer m, const integer n,
                           const dcomplex alpha, const dcomplex* a, const integer lda,
                                                 const dcomplex* x, const integer incx,
                           const dcomplex  beta,       dcomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(zgemv,ZGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_zgbmv(const char trans,
                           const integer m, const integer n, const integer kl, const integer ku,
                           const dcomplex alpha, const dcomplex* a, const integer lda,
                                                 const dcomplex* x, const integer incx,
                           const dcomplex  beta,       dcomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(zgbmv,ZGBMV)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_zhemv(const char uplo, const integer n,
                           const dcomplex alpha, const dcomplex* a, const integer lda,
                                                 const dcomplex* x, const integer incx,
                           const dcomplex  beta,       dcomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(zhemv,ZHEMV)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_zhbmv(const char uplo, const integer n, const integer k,
                           const dcomplex alpha, const dcomplex* a, const integer lda,
                                                 const dcomplex* x, const integer incx,
                           const dcomplex  beta,       dcomplex* y, const integer incy)
{
    MARRAY_FC_FUNC(zhbmv,ZHBMV)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_zhpmv(const char uplo, const integer n,
                           const dcomplex alpha, const dcomplex* ap,
                                                 const dcomplex*  x, const integer incx,
                           const dcomplex  beta,       dcomplex*  y, const integer incy)
{
    MARRAY_FC_FUNC(zhpmv,ZHPMV)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
}

static inline void c_ztrmv(const char uplo, const char trans, const char diag, const integer n,
                           const dcomplex* a, const integer lda,
                                 dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ztrmv,ZTRMV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_ztbmv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const dcomplex* a, const integer lda,
                                 dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ztbmv,ZTBMV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_ztpmv(const char uplo, const char trans, const char diag, const integer n,
                           const dcomplex* ap, dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ztpmv,ZTPMV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_ztrsv(const char uplo, const char trans, const char diag, const integer n,
                           const dcomplex* a, const integer lda,
                                 dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ztrsv,ZTRSV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_ztbsv(const char uplo, const char trans, const char diag,
                           const integer n, const integer k,
                           const dcomplex* a, const integer lda,
                                 dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ztbsv,ZTBSV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_ztpsv(const char uplo, const char trans, const char diag, const integer n,
                           const dcomplex* ap, dcomplex* x, const integer incx)
{
    MARRAY_FC_FUNC(ztpsv,ZTPSV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_zgeru(const integer m, const integer n,
                           const dcomplex alpha, const dcomplex* x, const integer incx,
                                                 const dcomplex* y, const integer incy,
                                                       dcomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(zgeru,ZGERU)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_zgerc(const integer m, const integer n,
                           const dcomplex alpha, const dcomplex* x, const integer incx,
                                                 const dcomplex* y, const integer incy,
                                                       dcomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(zgerc,ZGERC)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_zher (const char uplo, const integer n,
                           const double alpha, const dcomplex* x, const integer incx,
                                                       dcomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(zher,ZHER)(&uplo, &n, &alpha, x, &incx, a, &lda);
}

static inline void c_zhpr (const char uplo, const integer n,
                           const double alpha, const dcomplex* x, const integer incx,
                                                       dcomplex* ap)
{
    MARRAY_FC_FUNC(zhpr,ZHPR)(&uplo, &n, &alpha, x, &incx, ap);
}

static inline void c_zher2(const char uplo, const integer n,
                           const dcomplex alpha, const dcomplex* x, const integer incx,
                                                 const dcomplex* y, const integer incy,
                                                       dcomplex* a, const integer lda)
{
    MARRAY_FC_FUNC(zher2,ZHER2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_zhpr2(const char uplo, const integer n,
                           const dcomplex alpha, const dcomplex* x, const integer incx,
                                                 const dcomplex* y, const integer incy,
                                                       dcomplex* ap)
{
    MARRAY_FC_FUNC(zhpr2,ZHPR2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
}

/******************************************************************************
 *
 * Level 3 BLAS, C wrappers
 *
 *****************************************************************************/
static inline void c_sgemm (const char transa, const char transb,
                            const integer m, const integer n, const integer k,
                            const float alpha, const float* a, const integer lda,
                                               const float* b, const integer ldb,
                            const float  beta,       float* c, const integer ldc)
{
    MARRAY_FC_FUNC(sgemm,SGEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_ssymm (const char side, const char uplo,
                            const integer m, const integer n,
                            const float alpha, const float* a, const integer lda,
                                               const float* b, const integer ldb,
                            const float  beta,       float* c, const integer ldc)
{
    MARRAY_FC_FUNC(ssymm,SSYMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_ssyrk (const char uplo, const char trans,
                            const integer n, const integer k,
                            const float alpha, const float* a, const integer lda,
                            const float  beta,       float* c, const integer ldc)
{
    MARRAY_FC_FUNC(ssyrk,SSYRK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

static inline void c_ssyr2k(const char uplo, const char trans,
                            const integer n, const integer k,
                            const float alpha, const float* a, const integer lda,
                                               const float* b, const integer ldb,
                            const float  beta,       float* c, const integer ldc)
{
    MARRAY_FC_FUNC(ssyr2k,SSYR2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_strmm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const float alpha, const float* a, const integer lda,
                                                     float* b, const integer ldb)
{
    MARRAY_FC_FUNC(strmm,STRMM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_strsm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const float alpha, const float* a, const integer lda,
                                                     float* b, const integer ldb)
{
    MARRAY_FC_FUNC(strsm,STRSM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_dgemm (const char transa, const char transb,
                            const integer m, const integer n, const integer k,
                            const double alpha, const double* a, const integer lda,
                                                const double* b, const integer ldb,
                            const double  beta,       double* c, const integer ldc)
{
    MARRAY_FC_FUNC(dgemm,DGEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_dsymm (const char side, const char uplo,
                            const integer m, const integer n,
                            const double alpha, const double* a, const integer lda,
                                                const double* b, const integer ldb,
                            const double  beta,       double* c, const integer ldc)
{
    MARRAY_FC_FUNC(dsymm,DSYMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_dsyrk (const char uplo, const char trans,
                            const integer n, const integer k,
                            const double alpha, const double* a, const integer lda,
                            const double  beta,       double* c, const integer ldc)
{
    MARRAY_FC_FUNC(dsyrk,DSYRK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

static inline void c_dsyr2k(const char uplo, const char trans,
                            const integer n, const integer k,
                            const double alpha, const double* a, const integer lda,
                                                const double* b, const integer ldb,
                            const double  beta,       double* c, const integer ldc)
{
    MARRAY_FC_FUNC(dsyr2k,DSYR2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_dtrmm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const double alpha, const double* a, const integer lda,
                                                      double* b, const integer ldb)
{
    MARRAY_FC_FUNC(dtrmm,DTRMM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_dtrsm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const double alpha, const double* a, const integer lda,
                                                      double* b, const integer ldb)
{
    MARRAY_FC_FUNC(dtrsm,DTRSM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_cgemm (const char transa, const char transb,
                            const integer m, const integer n, const integer k,
                            const scomplex alpha, const scomplex* a, const integer lda,
                                                  const scomplex* b, const integer ldb,
                            const scomplex  beta,       scomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(cgemm,CGEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_chemm (const char side, const char uplo,
                            const integer m, const integer n,
                            const scomplex alpha, const scomplex* a, const integer lda,
                                                  const scomplex* b, const integer ldb,
                            const scomplex  beta,       scomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(chemm,CHEMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_csymm (const char side, const char uplo,
                            const integer m, const integer n,
                            const scomplex alpha, const scomplex* a, const integer lda,
                                                  const scomplex* b, const integer ldb,
                            const scomplex  beta,       scomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(csymm,CSYMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_csyrk (const char uplo, const char trans,
                            const integer n, const integer k,
                            const scomplex alpha, const scomplex* a, const integer lda,
                            const scomplex  beta,       scomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(csyrk,CSYRK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

static inline void c_csyr2k(const char uplo, const char trans,
                            const integer n, const integer k,
                            const scomplex alpha, const scomplex* a, const integer lda,
                                                  const scomplex* b, const integer ldb,
                            const scomplex  beta,       scomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(csyr2k,CSYR2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_cherk (const char uplo, const char trans,
                            const integer n, const integer k,
                            const float alpha, const scomplex* a, const integer lda,
                            const float  beta,       scomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(cherk,CHERK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

static inline void c_cher2k(const char uplo, const char trans,
                            const integer n, const integer k,
                            const scomplex alpha, const scomplex* a, const integer lda,
                                                  const scomplex* b, const integer ldb,
                            const    float  beta,       scomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(cher2k,CHER2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_ctrmm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const scomplex alpha, const scomplex* a, const integer lda,
                                                        scomplex* b, const integer ldb)
{
    MARRAY_FC_FUNC(ctrmm,CTRMM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_ctrsm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const scomplex alpha, const scomplex* a, const integer lda,
                                                        scomplex* b, const integer ldb)
{
    MARRAY_FC_FUNC(ctrsm,CTRSM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_zgemm (const char transa, const char transb,
                            const integer m, const integer n, const integer k,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                                                  const dcomplex* b, const integer ldb,
                            const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(zgemm,ZGEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_zhemm (const char side, const char uplo,
                            const integer m, const integer n,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                                                  const dcomplex* b, const integer ldb,
                            const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(zhemm,ZHEMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_zsymm (const char side, const char uplo,
                            const integer m, const integer n,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                                                  const dcomplex* b, const integer ldb,
                            const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(zsymm,ZSYMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_zsyrk (const char uplo, const char trans,
                            const integer n, const integer k,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                            const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(zsyrk,ZSYRK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

static inline void c_zsyr2k(const char uplo, const char trans,
                            const integer n, const integer k,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                                                  const dcomplex* b, const integer ldb,
                            const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(zsyr2k,ZSYR2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_zherk (const char uplo, const char trans,
                            const integer n, const integer k,
                            const double alpha, const dcomplex* a, const integer lda,
                            const double  beta,       dcomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(zherk,ZHERK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

static inline void c_zher2k(const char uplo, const char trans,
                            const integer n, const integer k,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                                                  const dcomplex* b, const integer ldb,
                            const   double  beta,       dcomplex* c, const integer ldc)
{
    MARRAY_FC_FUNC(zher2k,ZHER2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_ztrmm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                                                        dcomplex* b, const integer ldb)
{
    MARRAY_FC_FUNC(ztrmm,ZTRMM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_ztrsm (const char side, const char uplo, const char transa, const char diag,
                            const integer m, const integer n,
                            const dcomplex alpha, const dcomplex* a, const integer lda,
                                                        dcomplex* b, const integer ldb)
{
    MARRAY_FC_FUNC(ztrsm,ZTRSM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

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
namespace blas
{

namespace detail
{

using namespace MArray::detail;

template <typename T>
decltype(auto) view(T&& x)
{
    if constexpr (MArray::detail::is_marray_v<T> ||
                  MArray::detail::is_marray_slice_v<T>)
    {
        return x.view();
    }
    else
    {
        return std::forward<T>(x);
    }
}

#define MARRAY_FORWARD_AS_VIEW(name) \
template <typename... Args> \
auto name(Args&&... args) \
{ \
    name##_(detail::view(std::forward<Args>(args))...); \
}

}

/******************************************************************************
 *
 * Level 1 BLAS, C++ overloads
 *
 *****************************************************************************/
inline void rotg(float* a, float* b, float* c, float* s)
{
    srotg(a, b, c, s);
}

inline void rotg(double* a, double* b, double* c, double* s)
{
    drotg(a, b, c, s);
}

inline void rotg(scomplex* a, scomplex* b, float* c, scomplex* s)
{
    crotg(a, b, c, s);
}

inline void rotg(dcomplex* a, dcomplex* b, double* c, dcomplex* s)
{
    zrotg(a, b, c, s);
}

inline void rotmg(float* d1, float* d2, float* a, const float b, float* param)
{
    srotmg(d1, d2, a, b, param);
}

inline void rotmg(double* d1, double* d2, double* a, const double b, double* param)
{
    drotmg(d1, d2, a, b, param);
}

inline void rot(const integer n, float* x, const integer incx,
                                 float* y, const integer incy, const float c, const float s)
{
    srot(n, x, incx, y, incy, c, s);
}

inline void rot(const integer n, double* x, const integer incx,
                                 double* y, const integer incy, const double c, const double s)
{
    drot(n, x, incx, y, incy, c, s);
}

inline void rot(const integer n, scomplex* x, const integer incx,
                                 scomplex* y, const integer incy, const float c, const float s)
{
    csrot(n, x, incx, y, incy, c, s);
}

inline void rot(const integer n, dcomplex* x, const integer incx,
                                 dcomplex* y, const integer incy, const double c, const double s)
{
    zdrot(n, x, incx, y, incy, c, s);
}

inline void rotm(const integer n, float* x, const integer incx,
                                  float* y, const integer incy, float* param)
{
    srotm(n, x, incx, y, incy, param);
}

inline void rotm(const integer n, double* x, const integer incx,
                                  double* y, const integer incy, double* param)
{
    drotm(n, x, incx, y, incy, param);
}

inline void swap(const integer n, float* x, const integer incx,
                                  float* y, const integer incy)
{
    sswap(n, x, incx, y, incy);
}

inline void swap(const integer n, double* x, const integer incx,
                                  double* y, const integer incy)
{
    dswap(n, x, incx, y, incy);
}

inline void swap(const integer n, scomplex* x, const integer incx,
                                  scomplex* y, const integer incy)
{
    cswap(n, x, incx, y, incy);
}

inline void swap(const integer n, dcomplex* x, const integer incx,
                                  dcomplex* y, const integer incy)
{
    zswap(n, x, incx, y, incy);
}

inline void scal(const integer n, const float alpha, float* x, const integer incx)
{
    sscal(n, alpha, x, incx);
}

inline void scal(const integer n, const double alpha, double* x, const integer incx)
{
    dscal(n, alpha, x, incx);
}

inline void scal(const integer n, const scomplex alpha, scomplex* x, const integer incx)
{
    cscal(n, alpha, x, incx);
}

inline void scal(const integer n, const dcomplex alpha, dcomplex* x, const integer incx)
{
    zscal(n, alpha, x, incx);
}

inline void scal(const integer n, const float alpha, scomplex* x, const integer incx)
{
    csscal(n, alpha, x, incx);
}

inline void scal(const integer n, const double alpha, dcomplex* x, const integer incx)
{
    zdscal(n, alpha, x, incx);
}

inline void copy(const integer n, const float* x, const integer incx,
                                        float* y, const integer incy)
{
    scopy(n, x, incx, y, incy);
}

inline void copy(const integer n, const double* x, const integer incx,
                                        double* y, const integer incy)
{
    dcopy(n, x, incx, y, incy);
}

inline void copy(const integer n, const scomplex* x, const integer incx,
                                        scomplex* y, const integer incy)
{
    ccopy(n, x, incx, y, incy);
}

inline void copy(const integer n, const dcomplex* x, const integer incx,
                                        dcomplex* y, const integer incy)
{
    zcopy(n, x, incx, y, incy);
}

inline void axpy(const integer n, const float alpha, const float* x, const integer incx,
                                                           float* y, const integer incy)
{
    saxpy(n, alpha, x, incx, y, incy);
}

inline void axpy(const integer n, const double alpha, const double* x, const integer incx,
                                                            double* y, const integer incy)
{
    daxpy(n, alpha, x, incx, y, incy);
}

inline void axpy(const integer n, const scomplex alpha, const scomplex* x, const integer incx,
                                                              scomplex* y, const integer incy)
{
    caxpy(n, alpha, x, incx, y, incy);
}

inline void axpy(const integer n, const dcomplex alpha, const dcomplex* x, const integer incx,
                                                              dcomplex* y, const integer incy)
{
    zaxpy(n, alpha, x, incx, y, incy);
}

inline float dot(const integer n, const float* x, const integer incx,
                                  const float* y, const integer incy)
{
    return sdot(n, x, incx, y, incy);
}

inline double dot(const integer n, const double* x, const integer incx,
                                   const double* y, const integer incy)
{
    return ddot(n, x, incx, y, incy);
}

inline scomplex dot(const integer n, const scomplex* x, const integer incx,
                                     const scomplex* y, const integer incy)
{
    return cdotc(n, x, incx, y, incy);
}

inline dcomplex dot(const integer n, const dcomplex* x, const integer incx,
                                     const dcomplex* y, const integer incy)
{
    return zdotc(n, x, incx, y, incy);
}

inline float dotc(const integer n, const float* x, const integer incx,
                                   const float* y, const integer incy)
{
    return sdot(n, x, incx, y, incy);
}

inline double dotc(const integer n, const double* x, const integer incx,
                                    const double* y, const integer incy)
{
    return ddot(n, x, incx, y, incy);
}

inline scomplex dotc(const integer n, const scomplex* x, const integer incx,
                                      const scomplex* y, const integer incy)
{
    return cdotc(n, x, incx, y, incy);
}

inline dcomplex dotc(const integer n, const dcomplex* x, const integer incx,
                                      const dcomplex* y, const integer incy)
{
    return zdotc(n, x, incx, y, incy);
}

inline float dotu(const integer n, const float* x, const integer incx,
                                   const float* y, const integer incy)
{
    return sdot(n, x, incx, y, incy);
}

inline double dotu(const integer n, const double* x, const integer incx,
                                    const double* y, const integer incy)
{
    return ddot(n, x, incx, y, incy);
}

inline scomplex dotu(const integer n, const scomplex* x, const integer incx,
                                      const scomplex* y, const integer incy)
{
    return cdotu(n, x, incx, y, incy);
}

inline dcomplex dotu(const integer n, const dcomplex* x, const integer incx,
                                      const dcomplex* y, const integer incy)
{
    return zdotu(n, x, incx, y, incy);
}

inline float nrm2(const integer n, const float* x, const integer incx)
{
    return snrm2(n, x, incx);
}

inline double nrm2(const integer n, const double* x, const integer incx)
{
    return dnrm2(n, x, incx);
}

inline float nrm2(const integer n, const scomplex* x, const integer incx)
{
    return scnrm2(n, x, incx);
}

inline double nrm2(const integer n, const dcomplex* x, const integer incx)
{
    return dznrm2(n, x, incx);
}

inline float asum(const integer n, const float* x, const integer incx)
{
    return sasum(n, x, incx);
}

inline double asum(const integer n, const double* x, const integer incx)
{
    return dasum(n, x, incx);
}

inline float asum(const integer n, const scomplex* x, const integer incx)
{
    return scasum(n, x, incx);
}

inline double asum(const integer n, const dcomplex* x, const integer incx)
{
    return dzasum(n, x, incx);
}

inline integer iamax(const integer n, const float* x, const integer incx)
{
    return isamax(n, x, incx);
}

inline integer iamax(const integer n, const double* x, const integer incx)
{
    return idamax(n, x, incx);
}

inline integer iamax(const integer n, const scomplex* x, const integer incx)
{
    return icamax(n, x, incx);
}

inline integer iamax(const integer n, const dcomplex* x, const integer incx)
{
    return izamax(n, x, incx);
}

inline float amax(const integer n, const float* x, const integer incx)
{
    return std::abs(x[isamax(n, x, incx)]);
}

inline double amax(const integer n, const double* x, const integer incx)
{
    return std::abs(x[idamax(n, x, incx)]);
}

inline float amax(const integer n, const scomplex* x, const integer incx)
{
    return std::abs(x[icamax(n, x, incx)]);
}

inline double amax(const integer n, const dcomplex* x, const integer incx)
{
    return std::abs(x[izamax(n, x, incx)]);
}

/******************************************************************************
 *
 * Level 2 BLAS, C++ overloads
 *
 *****************************************************************************/
inline void gemv(const char trans, const integer m, const integer n,
                 const float alpha, const float* a, const integer lda,
                                    const float* x, const integer incx,
                 const float  beta,       float* y, const integer incy)
{
    sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void gemv(const char trans, const integer m, const integer n,
                 const double alpha, const double* a, const integer lda,
                                     const double* x, const integer incx,
                 const double  beta,       double* y, const integer incy)
{
    dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void gemv(const char trans, const integer m, const integer n,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                       const scomplex* x, const integer incx,
                 const scomplex  beta,       scomplex* y, const integer incy)
{
    cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void gemv(const char trans, const integer m, const integer n,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                       const dcomplex* x, const integer incx,
                 const dcomplex  beta,       dcomplex* y, const integer incy)
{
    zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void gbmv(const char trans,
                 const integer m, const integer n, const integer kl, const integer ku,
                 const float alpha, const float* a, const integer lda,
                                    const float* x, const integer incx,
                 const float  beta,       float* y, const integer incy)
{
    sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}

inline void gbmv(const char trans,
                 const integer m, const integer n, const integer kl, const integer ku,
                 const double alpha, const double* a, const integer lda,
                                     const double* x, const integer incx,
                 const double  beta,       double* y, const integer incy)
{
    dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}

inline void gbmv(const char trans,
                 const integer m, const integer n, const integer kl, const integer ku,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                       const scomplex* x, const integer incx,
                 const scomplex  beta,       scomplex* y, const integer incy)
{
    cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}

inline void gbmv(const char trans,
                 const integer m, const integer n, const integer kl, const integer ku,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                       const dcomplex* x, const integer incx,
                 const dcomplex  beta,       dcomplex* y, const integer incy)
{
    zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hemv(const char uplo, const integer n,
                 const float alpha, const float* a, const integer lda,
                                    const float* x, const integer incx,
                 const float  beta,       float* y, const integer incy)
{
    ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hemv(const char uplo, const integer n,
                 const double alpha, const double* a, const integer lda,
                                     const double* x, const integer incx,
                 const double  beta,       double* y, const integer incy)
{
    dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hemv(const char uplo, const integer n,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                       const scomplex* x, const integer incx,
                 const scomplex  beta,       scomplex* y, const integer incy)
{
    chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hemv(const char uplo, const integer n,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                       const dcomplex* x, const integer incx,
                 const dcomplex  beta,       dcomplex* y, const integer incy)
{
    zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hbmv(const char uplo, const integer n, const integer k,
                 const float alpha, const float* a, const integer lda,
                                    const float* x, const integer incx,
                 const float  beta,       float* y, const integer incy)
{
    ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hbmv(const char uplo, const integer n, const integer k,
                 const double alpha, const double* a, const integer lda,
                                     const double* x, const integer incx,
                 const double  beta,       double* y, const integer incy)
{
    dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hbmv(const char uplo, const integer n, const integer k,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                       const scomplex* x, const integer incx,
                 const scomplex  beta,       scomplex* y, const integer incy)
{
    chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hbmv(const char uplo, const integer n, const integer k,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                       const dcomplex* x, const integer incx,
                 const dcomplex  beta,       dcomplex* y, const integer incy)
{
    zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
}

inline void hpmv(const char uplo, const integer n,
                 const float alpha, const float* ap,
                                    const float*  x, const integer incx,
                 const float  beta,       float*  y, const integer incy)
{
    sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy);
}

inline void hpmv(const char uplo, const integer n,
                 const double alpha, const double* ap,
                                     const double*  x, const integer incx,
                 const double  beta,       double*  y, const integer incy)
{
    dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy);
}

inline void hpmv(const char uplo, const integer n,
                 const scomplex alpha, const scomplex* ap,
                                       const scomplex*  x, const integer incx,
                 const scomplex  beta,       scomplex*  y, const integer incy)
{
    chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy);
}

inline void hpmv(const char uplo, const integer n,
                 const dcomplex alpha, const dcomplex* ap,
                                       const dcomplex*  x, const integer incx,
                 const dcomplex  beta,       dcomplex*  y, const integer incy)
{
    zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy);
}

inline void trmv(const char uplo, const char trans, const char diag, const integer n,
                 const float* a, const integer lda,
                       float* x, const integer incx)
{
    strmv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void trmv(const char uplo, const char trans, const char diag, const integer n,
                 const double* a, const integer lda,
                       double* x, const integer incx)
{
    dtrmv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void trmv(const char uplo, const char trans, const char diag, const integer n,
                 const scomplex* a, const integer lda,
                       scomplex* x, const integer incx)
{
    ctrmv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void trmv(const char uplo, const char trans, const char diag, const integer n,
                 const dcomplex* a, const integer lda,
                       dcomplex* x, const integer incx)
{
    ztrmv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void tbmv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const float* a, const integer lda,
                       float* x, const integer incx)
{
    stbmv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tbmv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const double* a, const integer lda,
                       double* x, const integer incx)
{
    dtbmv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tbmv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const scomplex* a, const integer lda,
                       scomplex* x, const integer incx)
{
    ctbmv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tbmv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const dcomplex* a, const integer lda,
                       dcomplex* x, const integer incx)
{
    ztbmv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tpmv(const char uplo, const char trans, const char diag, const integer n,
                 const float* ap, float* x, const integer incx)
{
    stpmv(uplo, trans, diag, n, ap, x, incx);
}

inline void tpmv(const char uplo, const char trans, const char diag, const integer n,
                 const double* ap, double* x, const integer incx)
{
    dtpmv(uplo, trans, diag, n, ap, x, incx);
}

inline void tpmv(const char uplo, const char trans, const char diag, const integer n,
                 const scomplex* ap, scomplex* x, const integer incx)
{
    ctpmv(uplo, trans, diag, n, ap, x, incx);
}

inline void tpmv(const char uplo, const char trans, const char diag, const integer n,
                 const dcomplex* ap, dcomplex* x, const integer incx)
{
    ztpmv(uplo, trans, diag, n, ap, x, incx);
}

inline void trsv(const char uplo, const char trans, const char diag, const integer n,
                 const float* a, const integer lda,
                       float* x, const integer incx)
{
    strsv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void trsv(const char uplo, const char trans, const char diag, const integer n,
                 const double* a, const integer lda,
                       double* x, const integer incx)
{
    dtrsv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void trsv(const char uplo, const char trans, const char diag, const integer n,
                 const scomplex* a, const integer lda,
                       scomplex* x, const integer incx)
{
    ctrsv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void trsv(const char uplo, const char trans, const char diag, const integer n,
                 const dcomplex* a, const integer lda,
                       dcomplex* x, const integer incx)
{
    ztrsv(uplo, trans, diag, n, a, lda, x, incx);
}

inline void tbsv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const float* a, const integer lda,
                       float* x, const integer incx)
{
    stbsv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tbsv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const double* a, const integer lda,
                       double* x, const integer incx)
{
    dtbsv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tbsv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const scomplex* a, const integer lda,
                       scomplex* x, const integer incx)
{
    ctbsv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tbsv(const char uplo, const char trans, const char diag,
                 const integer n, const integer k,
                 const dcomplex* a, const integer lda,
                       dcomplex* x, const integer incx)
{
    ztbsv(uplo, trans, diag, n, k, a, lda, x, incx);
}

inline void tpsv(const char uplo, const char trans, const char diag, const integer n,
                 const float* ap, float* x, const integer incx)
{
    stpsv(uplo, trans, diag, n, ap, x, incx);
}

inline void tpsv(const char uplo, const char trans, const char diag, const integer n,
                 const double* ap, double* x, const integer incx)
{
    dtpsv(uplo, trans, diag, n, ap, x, incx);
}

inline void tpsv(const char uplo, const char trans, const char diag, const integer n,
                 const scomplex* ap, scomplex* x, const integer incx)
{
    ctpsv(uplo, trans, diag, n, ap, x, incx);
}

inline void tpsv(const char uplo, const char trans, const char diag, const integer n,
                const dcomplex* ap, dcomplex* x, const integer incx)
{
    ztpsv(uplo, trans, diag, n, ap, x, incx);
}

inline void ger(const integer m, const integer n,
                const float alpha, const float* x, const integer incx,
                                   const float* y, const integer incy,
                                         float* a, const integer lda)
{
    sger(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void ger(const integer m, const integer n,
                const double alpha, const double* x, const integer incx,
                                    const double* y, const integer incy,
                                          double* a, const integer lda)
{
    dger(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void ger(const integer m, const integer n,
                const scomplex alpha, const scomplex* x, const integer incx,
                                      const scomplex* y, const integer incy,
                                            scomplex* a, const integer lda)
{
    cgerc(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void ger(const integer m, const integer n,
                const dcomplex alpha, const dcomplex* x, const integer incx,
                                      const dcomplex* y, const integer incy,
                                            dcomplex* a, const integer lda)
{
    zgerc(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void gerc(const integer m, const integer n,
                 const float alpha, const float* x, const integer incx,
                                    const float* y, const integer incy,
                                          float* a, const integer lda)
{
    sger(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void gerc(const integer m, const integer n,
                 const double alpha, const double* x, const integer incx,
                                     const double* y, const integer incy,
                                           double* a, const integer lda)
{
    dger(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void gerc(const integer m, const integer n,
                 const scomplex alpha, const scomplex* x, const integer incx,
                                       const scomplex* y, const integer incy,
                                             scomplex* a, const integer lda)
{
    cgerc(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void gerc(const integer m, const integer n,
                 const dcomplex alpha, const dcomplex* x, const integer incx,
                                       const dcomplex* y, const integer incy,
                                             dcomplex* a, const integer lda)
{
    zgerc(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void geru(const integer m, const integer n,
                 const float alpha, const float* x, const integer incx,
                                    const float* y, const integer incy,
                                          float* a, const integer lda)
{
    sger(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void geru(const integer m, const integer n,
                 const double alpha, const double* x, const integer incx,
                                     const double* y, const integer incy,
                                           double* a, const integer lda)
{
    dger(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void geru(const integer m, const integer n,
                 const scomplex alpha, const scomplex* x, const integer incx,
                                       const scomplex* y, const integer incy,
                                             scomplex* a, const integer lda)
{
    cgeru(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void geru(const integer m, const integer n,
                 const dcomplex alpha, const dcomplex* x, const integer incx,
                                       const dcomplex* y, const integer incy,
                                             dcomplex* a, const integer lda)
{
    zgeru(m, n, alpha, x, incx, y, incy, a, lda);
}

inline void her(const char uplo, const integer n,
                const float alpha, const float* x, const integer incx,
                                         float* a, const integer lda)
{
    ssyr(uplo, n, alpha, x, incx, a, lda);
}

inline void her(const char uplo, const integer n,
                const double alpha, const double* x, const integer incx,
                                          double* a, const integer lda)
{
    dsyr(uplo, n, alpha, x, incx, a, lda);
}

inline void her(const char uplo, const integer n,
                const float alpha, const scomplex* x, const integer incx,
                                            scomplex* a, const integer lda)
{
    cher(uplo, n, alpha, x, incx, a, lda);
}

inline void her(const char uplo, const integer n,
                const double alpha, const dcomplex* x, const integer incx,
                                            dcomplex* a, const integer lda)
{
    zher(uplo, n, alpha, x, incx, a, lda);
}

inline void hpr(const char uplo, const integer n,
                const float alpha, const float* x, const integer incx,
                                         float* ap)
{
    sspr(uplo, n, alpha, x, incx, ap);
}

inline void hpr(const char uplo, const integer n,
                const double alpha, const double* x, const integer incx,
                                          double* ap)
{
    dspr(uplo, n, alpha, x, incx, ap);
}

inline void hpr(const char uplo, const integer n,
                const float alpha, const scomplex* x, const integer incx,
                                            scomplex* ap)
{
    chpr(uplo, n, alpha, x, incx, ap);
}

inline void hpr(const char uplo, const integer n,
                const double alpha, const dcomplex* x, const integer incx,
                                            dcomplex* ap)
{
    zhpr(uplo, n, alpha, x, incx, ap);
}

inline void her2(const char uplo, const integer n,
                 const float alpha, const float* x, const integer incx,
                                    const float* y, const integer incy,
                                          float* a, const integer lda)
{
    ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda);
}

inline void her2(const char uplo, const integer n,
                 const double alpha, const double* x, const integer incx,
                                     const double* y, const integer incy,
                                           double* a, const integer lda)
{
    dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda);
}
inline void her2(const char uplo, const integer n,
                 const scomplex alpha, const scomplex* x, const integer incx,
                                       const scomplex* y, const integer incy,
                                             scomplex* a, const integer lda)
{
    cher2(uplo, n, alpha, x, incx, y, incy, a, lda);
}

inline void her2(const char uplo, const integer n,
                 const dcomplex alpha, const dcomplex* x, const integer incx,
                                       const dcomplex* y, const integer incy,
                                             dcomplex* a, const integer lda)
{
    zher2(uplo, n, alpha, x, incx, y, incy, a, lda);
}

inline void hpr2(const char uplo, const integer n,
                 const float alpha, const float* x, const integer incx,
                                    const float* y, const integer incy,
                                          float* ap)
{
    sspr2(uplo, n, alpha, x, incx, y, incy, ap);
}

inline void hpr2(const char uplo, const integer n,
                 const double alpha, const double* x, const integer incx,
                                     const double* y, const integer incy,
                                           double* ap)
{
    dspr2(uplo, n, alpha, x, incx, y, incy, ap);
}

inline void hpr2(const char uplo, const integer n,
                 const scomplex alpha, const scomplex* x, const integer incx,
                                       const scomplex* y, const integer incy,
                                             scomplex* ap)
{
    chpr2(uplo, n, alpha, x, incx, y, incy, ap);
}

inline void hpr2(const char uplo, const integer n,
                 const dcomplex alpha, const dcomplex* x, const integer incx,
                                       const dcomplex* y, const integer incy,
                                             dcomplex* ap)
{
    zhpr2(uplo, n, alpha, x, incx, y, incy, ap);
}

/******************************************************************************
 *
 * Level 3 BLAS, C++ overloads
 *
 *****************************************************************************/
inline void gemm(const char transa, const char transb,
                 const integer m, const integer n, const integer k,
                 const float alpha, const float* a, const integer lda,
                                    const float* b, const integer ldb,
                 const float  beta,       float* c, const integer ldc)
{
    sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void gemm(const char transa, const char transb,
                 const integer m, const integer n, const integer k,
                 const double alpha, const double* a, const integer lda,
                                     const double* b, const integer ldb,
                 const double  beta,       double* c, const integer ldc)
{
    dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void gemm(const char transa, const char transb,
                 const integer m, const integer n, const integer k,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                       const scomplex* b, const integer ldb,
                 const scomplex  beta,       scomplex* c, const integer ldc)
{
    cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void gemm(const char transa, const char transb,
                 const integer m, const integer n, const integer k,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                       const dcomplex* b, const integer ldb,
                 const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void hemm(const char side, const char uplo,
                 const integer m, const integer n,
                 const float alpha, const float* a, const integer lda,
                                    const float* b, const integer ldb,
                 const float  beta,       float* c, const integer ldc)
{
    ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void hemm(const char side, const char uplo,
                 const integer m, const integer n,
                 const double alpha, const double* a, const integer lda,
                                     const double* b, const integer ldb,
                 const double  beta,       double* c, const integer ldc)
{
    dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void hemm(const char side, const char uplo,
                 const integer m, const integer n,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                       const scomplex* b, const integer ldb,
                 const scomplex  beta,       scomplex* c, const integer ldc)
{
    chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void hemm(const char side, const char uplo,
                 const integer m, const integer n,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                       const dcomplex* b, const integer ldb,
                 const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void symm(const char side, const char uplo,
                 const integer m, const integer n,
                 const float alpha, const float* a, const integer lda,
                                    const float* b, const integer ldb,
                 const float  beta,       float* c, const integer ldc)
{
    ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void symm(const char side, const char uplo,
                 const integer m, const integer n,
                 const double alpha, const double* a, const integer lda,
                                     const double* b, const integer ldb,
                 const double  beta,       double* c, const integer ldc)
{
    dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void symm(const char side, const char uplo,
                 const integer m, const integer n,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                       const scomplex* b, const integer ldb,
                 const scomplex  beta,       scomplex* c, const integer ldc)
{
    csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void symm(const char side, const char uplo,
                 const integer m, const integer n,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                       const dcomplex* b, const integer ldb,
                 const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void syrk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const float alpha, const float* a, const integer lda,
                 const float  beta,       float* c, const integer ldc)
{
    ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void syrk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const double alpha, const double* a, const integer lda,
                 const double  beta,       double* c, const integer ldc)
{
    dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void syrk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const scomplex alpha, const scomplex* a, const integer lda,
                 const scomplex  beta,       scomplex* c, const integer ldc)
{
    csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void syrk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                 const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void herk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const float alpha, const float* a, const integer lda,
                 const float  beta,       float* c, const integer ldc)
{
    ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void herk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const double alpha, const double* a, const integer lda,
                 const double  beta,       double* c, const integer ldc)
{
    dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void herk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const float alpha, const scomplex* a, const integer lda,
                 const float  beta,       scomplex* c, const integer ldc)
{
    cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void herk(const char uplo, const char trans,
                 const integer n, const integer k,
                 const double alpha, const dcomplex* a, const integer lda,
                 const double  beta,       dcomplex* c, const integer ldc)
{
    zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

inline void syr2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const float alpha, const float* a, const integer lda,
                                     const float* b, const integer ldb,
                  const float  beta,       float* c, const integer ldc)
{
    ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void syr2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const double alpha, const double* a, const integer lda,
                                      const double* b, const integer ldb,
                  const double  beta,       double* c, const integer ldc)
{
    dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void syr2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const scomplex alpha, const scomplex* a, const integer lda,
                                        const scomplex* b, const integer ldb,
                  const scomplex  beta,       scomplex* c, const integer ldc)
{
    csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void syr2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const dcomplex alpha, const dcomplex* a, const integer lda,
                                        const dcomplex* b, const integer ldb,
                  const dcomplex  beta,       dcomplex* c, const integer ldc)
{
    zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void her2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const float alpha, const float* a, const integer lda,
                                     const float* b, const integer ldb,
                  const float  beta,       float* c, const integer ldc)
{
    ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void her2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const double alpha, const double* a, const integer lda,
                                      const double* b, const integer ldb,
                  const double  beta,       double* c, const integer ldc)
{
    dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void her2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const scomplex alpha, const scomplex* a, const integer lda,
                                        const scomplex* b, const integer ldb,
                  const    float  beta,       scomplex* c, const integer ldc)
{
    cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void her2k(const char uplo, const char trans,
                  const integer n, const integer k,
                  const dcomplex alpha, const dcomplex* a, const integer lda,
                                        const dcomplex* b, const integer ldb,
                  const   double  beta,       dcomplex* c, const integer ldc)
{
    zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void trmm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const float alpha, const float* a, const integer lda,
                                          float* b, const integer ldb)
{
    strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

inline void trmm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const double alpha, const double* a, const integer lda,
                                           double* b, const integer ldb)
{
    dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

inline void trmm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                             scomplex* b, const integer ldb)
{
    ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

inline void trmm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                             dcomplex* b, const integer ldb)
{
    ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

inline void trsm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const float alpha, const float* a, const integer lda,
                                          float* b, const integer ldb)
{
    strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

inline void trsm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const double alpha, const double* a, const integer lda,
                                           double* b, const integer ldb)
{
    dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

inline void trsm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const scomplex alpha, const scomplex* a, const integer lda,
                                             scomplex* b, const integer ldb)
{
    ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

inline void trsm(const char side, const char uplo, const char transa, const char diag,
                 const integer m, const integer n,
                 const dcomplex alpha, const dcomplex* a, const integer lda,
                                             dcomplex* b, const integer ldb)
{
    ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

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
    swap(x.length(), x.data(), x.stride(), y.data(), y.stride());
}

template <typename T>
std::enable_if_t<detail::is_marray_like_v<T,1>>
conj(T&& x_)
{
    auto x = x_.view();
    auto n = x.length();
    auto ptr = x.data();
    auto stride = x.stride();
    for (len_type i = 0;i < n;i++)
        ptr[i*stride] = std::conj(ptr[i*stride]);
}

template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<U,1>>
scal(T alpha, U&& x_)
{
    auto x = x_.view();
    scal(x.length(), alpha, x.data(), x.stride());
}

template <typename T, typename U>
std::enable_if_t<detail::is_marray_like_v<T,1> &&
                 detail::is_marray_like_v<U,1>>
copy(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());
    copy(x.length(), x.data(), x.stride(), y.data(), y.stride());
}

template <typename T, typename U, typename V>
std::enable_if_t<detail::is_marray_like_v<U,1> &&
                 detail::is_marray_like_v<V,1>>
axpy(T alpha, U&& x_, V&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());
    axpy(x.length(), alpha, x.data(), x.stride(), y.data(), y.stride());
}

template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto dotu(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());
    return dotu(x.length(), x.data(), x.stride(), y.data(), y.stride());
}

template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto dotc(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());
    return dotc(x.length(), x.data(), x.stride(), y.data(), y.stride());
}

template <typename T, typename U, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1> &&
                     detail::is_marray_like_v<U,1>>>
auto dot(T&& x_, U&& y_)
{
    auto x = x_.view();
    auto y = y_.view();
    MARRAY_ASSERT(x.length() == y.length());
    return dot(x.length(), x.data(), x.stride(), y.data(), y.stride());
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto nrm2(T&& x_)
{
    auto x = x_.view();
    return nrm2(x.length(), x.data(), x.stride());
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto asum(T&& x_)
{
    auto x = x_.view();
    return asum(x.length(), x.data(), x.stride());
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto iamax(T&& x_)
{
    auto x = x_.view();
    return iamax(x.length(), x.data(), x.stride());
}

template <typename T, typename=
    std::enable_if_t<detail::is_marray_like_v<T,1>>>
auto amax(T&& x_)
{
    auto x = x_.view();
    return amax(x.length(), x.data(), x.stride());
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
    marray<std::decay_t<typename U::value_type>,1> y({A.length(0)}, uninitialized);
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

    if (A_.stride(0) > 1)
    {
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;

        if constexpr (detail::is_complex_v<std::decay_t<typename V::value_type>>)
        {
            // Compute (alpha Ax + beta y)^H = alpha' x^H A + beta' y^H

            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<std::decay_t<typename V::value_type>,1> conjx(x_);

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
    marray<std::decay_t<typename U::value_type>,1> y({A.length(0)}, uninitialized);
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
                x_.data(), x_.stride(),
          beta, y_.data(), y_.stride());
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
    marray<std::decay_t<typename U::value_type>,1> y({A.length(0)}, uninitialized);
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

    if (A_.stride(0) > 1 && 0)
    {
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;

        if constexpr (detail::is_complex_v<std::decay_t<typename U::value_type>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<std::decay_t<typename U::value_type>,1> conjx(x_);
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
    marray<std::decay_t<typename U::value_type>,2> A({x.length(), x.length()}, uninitialized);
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

    if (A_.stride(0) > 1)
    {
        A_.transpose();
        uplo = uplo == 'L' ? 'U' :
               uplo == 'U' ? 'L' : uplo;

        if constexpr (detail::is_complex_v<std::decay_t<typename U::value_type>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<std::decay_t<typename U::value_type>,1> conjx(x_);
            marray<std::decay_t<typename V::value_type>,1> conjy(y_);
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
    marray<std::decay_t<typename U::value_type>,2> A({x.length(), x.length()}, uninitialized);
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

    if (A_.stride(0) > 1)
    {
        A_.transpose();
        x_.swap(y_);

        if constexpr (detail::is_complex_v<std::decay_t<typename V::value_type>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<std::decay_t<typename V::value_type>,1> conjx(x_);
            conj(conjx);

            geru(A_.length(0), A_.length(1),
                 alpha, x_.data(), x_.stride(),
                        y_.data(), y_.stride(),
                        A_.data(), A_.stride(1));

            return;
        }
    }

    MARRAY_ASSERT(A_.stride(0) == 1);

    ger(A_.length(0), A_.length(1),
        alpha, x_.data(), x_.stride(),
               y_.data(), y_.stride(),
               A_.data(), A_.stride(1));
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
    marray<std::decay_t<typename U::value_type>,2> A({x.length(), y.length()}, uninitialized);
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

    if (A_.stride(0) > 1)
    {
        A_.transpose();
        x_.swap(y_);

        if constexpr (detail::is_complex_v<std::decay_t<typename V::value_type>>)
        {
            MARRAY_ASSERT(A_.stride(0) == 1);

            marray<std::decay_t<typename V::value_type>,1> conjx(x_);
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
    marray<std::decay_t<typename U::value_type>,2> A({x.length(), y.length()}, uninitialized);
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
    marray<std::decay_t<typename U::value_type>,2> A({x.length(), y.length()}, uninitialized);
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
                 detail::is_marray_like_v<W,2>>
gemm(T alpha, U&& A, V&& B, W beta, X&& C)
{
    auto A_ = A.view();
    auto B_ = B.view();
    auto C_ = C.view();

    auto m = C_.length(0);
    auto n = C_.length(1);
    auto k = A_.length(1);

    ALWAYS_ASSERT(A_.length(0) == m);
    ALWAYS_ASSERT(B_.length(1) == n);
    ALWAYS_ASSERT(B_.length(0) == k);

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
    marray<std::decay_t<typename U::value_type>,2> C({A.length(0), B.length(1)}, uninitialized);
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

    ALWAYS_ASSERT(B_.length(0) == m);
    ALWAYS_ASSERT(B_.length(1) == n);
    ALWAYS_ASSERT(A_.length(0) == side == 'L' ? m : n);
    ALWAYS_ASSERT(A_.length(1) == side == 'L' ? m : n);

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
    marray<std::decay_t<typename U::value_type>,2> C({B.length(0), B.length(1)}, uninitialized);
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

    ALWAYS_ASSERT(B_.length(0) == m);
    ALWAYS_ASSERT(B_.length(1) == n);
    ALWAYS_ASSERT(A_.length(0) == side == 'L' ? m : n);
    ALWAYS_ASSERT(A_.length(1) == side == 'L' ? m : n);

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
    marray<std::decay_t<typename U::value_type>,2> C({B.length(0), B.length(1)}, uninitialized);
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

    ALWAYS_ASSERT(A_.length(0) == side == 'L' ? m : n);
    ALWAYS_ASSERT(A_.length(1) == side == 'L' ? m : n);

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

    ALWAYS_ASSERT(A_.length(0) == side == 'L' ? m : n);
    ALWAYS_ASSERT(A_.length(1) == side == 'L' ? m : n);

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
    marray<std::decay_t<typename U::value_type>,2> C({A.length(0), A.length(0)}, uninitialized);
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
    marray<std::decay_t<typename U::value_type>,2> C({A.length(0), B.length(0)}, uninitialized);
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
    marray<std::decay_t<typename U::value_type>,2> C({A.length(0), A.length(0)}, uninitialized);
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
    marray<std::decay_t<typename U::value_type>,2> C({A.length(0), B.length(0)}, uninitialized);
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

} //namespace blas
} //namespace MArray

#endif //__cplusplus

#endif //MARRAY_BLAS_H
