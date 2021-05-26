/*
  ARPACK++ v1.2 2/20/2000
  c++ interface to ARPACK code.

  MODULE blas1f.h
  BLAS 1 and BLAS 2 FORTRAN routines.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef BLAS1F_H
#define BLAS1F_H

#include "arch.h"

extern "C"
{

  // Single precision real routines.

  float F77NAME(sasum)(const ARint *n, const float *dx, const ARint *incx);

  void F77NAME(saxpy)(const ARint *n, const float *da, const float *dx,
                      const ARint *incx, float *dy, const ARint *incy);

  void F77NAME(scopy)(const ARint *n, const float *dx, const ARint *incx,
                      float *dy, const ARint *incy);

  float F77NAME(sdot)(const ARint *n, const float *dx, const ARint *incx,
                      const float *dy, const ARint *incy);

  float F77NAME(snrm2)(const ARint *n, const float *dx, const ARint *incx);

  void F77NAME(srot)(const ARint *n, float *dx, const ARint *incx, float *dy,
                     const ARint *incy, const float *c, const float *s);

  void F77NAME(srotg)(float *da, float *db, float *c, float *s);

  void F77NAME(sscal)(const ARint *n, float *da, float *dx, const ARint *incx);

  void F77NAME(sswap)(const ARint *n, float *dx, const ARint *incx,
                      float *dy, const ARint *incy);

  ARint F77NAME(isamax)(const ARint *n, const float *dx, const ARint *incx);

  void F77NAME(sgemv)(const char* trans, const ARint *m, const ARint *n, 
                      const float *alpha, const float *a, const ARint *lda, 
                      const float *x, const ARint *incx, const float *beta, 
                      float *y, const ARint *incy);

  void F77NAME(sgbmv)(const char* trans, const ARint *m, const ARint *n, 
                      const ARint *kl, const ARint *ku, const float *alpha,
                      const float *a, const ARint *lda, const float *x,
                      const ARint *incx, const float *beta, float *y,
                      const ARint *incy);

  void F77NAME(ssbmv)(const char* uplo, const ARint *n, const ARint *k, 
                      const float *alpha, const float *a, const ARint *lda, 
                      const float *x, const ARint *incx, const float *beta, 
                      float *y, const ARint *incy);

// Double precision real routines.

  double F77NAME(dasum)(const ARint *n, const double *dx, const ARint *incx);

  void F77NAME(daxpy)(const ARint *n, const double *da, const double *dx,
                      const ARint *incx, double *dy, const ARint *incy);

  void F77NAME(dcopy)(const ARint *n, const double *dx, const ARint *incx,
                      double *dy, const ARint *incy);

  double F77NAME(ddot)(const ARint *n, const double *dx, const ARint *incx,
                       const double *dy, const ARint *incy);

  double F77NAME(dnrm2)(const ARint *n, const double *dx, const ARint *incx);

  void F77NAME(drot)(const ARint *n, double *dx, const ARint *incx, double *dy,
                     const ARint *incy, const double *c, const double *s);

  void F77NAME(drotg)(double *da, double *db, double *c, double *s);

  void F77NAME(dscal)(const ARint *n, double *da, double *dx, const ARint *incx);

  void F77NAME(dswap)(const ARint *n, double *dx, const ARint *incx,
                      double *dy, const ARint *incy);

  ARint F77NAME(idamax)(const ARint *n, const double *dx, const ARint *incx);

  void F77NAME(dgemv)(const char* trans, const ARint *m, const ARint *n, 
                      const double *alpha, const double *a, const ARint *lda,
                      const double *x, const ARint *incx, const double *beta,
                      double *y, const ARint *incy);

  void F77NAME(dgbmv)(const char* trans, const ARint *m, const ARint *n, 
                      const ARint *kl, const ARint *ku, const double *alpha,
                      const double *a, const ARint *lda, const double *x,
                      const ARint *incx, const double *beta, double *y,
                      const ARint *incy);

  void F77NAME(dsbmv)(const char* uplo, const ARint *n, const ARint *k, 
                      const double *alpha, const double *a, const ARint *lda, 
                      const double *x, const ARint *incx, const double *beta, 
                      double *y, const ARint *incy);

  // Single precision complex routines.

#ifdef ARCOMP_H

  void F77NAME(cdotc)(arcomplex<float> *c, const ARint *n,
                      const arcomplex<float> *cx, const ARint *incx,
                      const arcomplex<float> *cy, const ARint *incy);

  void F77NAME(cdotu)(arcomplex<float> *c, const ARint *n,
                      const arcomplex<float> *cx, const ARint *incx,
                      const arcomplex<float> *cy, const ARint *incy);

  void F77NAME(caxpy)(const ARint *n, const arcomplex<float> *da,
                      const arcomplex<float> *dx, const ARint *incx,
                      arcomplex<float> *dy, const ARint *incy);

  void F77NAME(ccopy)(const ARint *n, const arcomplex<float> *dx,
                      const ARint *incx, arcomplex<float> *dy,
                      const ARint *incy);

  float F77NAME(scasum)(const ARint *n, const arcomplex<float> *dx,
                        const ARint *incx);

  float F77NAME(scnrm2)(const ARint *n, const arcomplex<float> *dx,
                        const ARint *incx);

  void F77NAME(csscal)(const ARint *n, const float *da, arcomplex<float> *dx,
                       const ARint *incx);

  void F77NAME(cscal)(const ARint *n, const arcomplex<float> *da,
                      arcomplex<float> *dx, const ARint *incx);

  ARint F77NAME(icamax)(const ARint *n, const arcomplex<float> *dx,
                          const ARint *incx);

  void F77NAME(cswap)(const ARint *n, arcomplex<float> *dx, 
                      const ARint *incx, arcomplex<float> *dy, 
                      const ARint *incy);

  void F77NAME(cgemv)(const char* trans, const ARint *m, 
                      const ARint *n, const arcomplex<float> *alpha,
                      const arcomplex<float> *a, const ARint *lda, 
                      const arcomplex<float> *x, const ARint *incx, 
                      const arcomplex<float> *beta, arcomplex<float> *y,
                      const ARint *incy);

  void F77NAME(cgbmv)(const char* trans, const ARint *m, 
                      const ARint *n, const ARint *kl, 
                      const ARint *ku, const arcomplex<float> *alpha,
                      const arcomplex<float> *a, const ARint *lda, 
                      const arcomplex<float> *x, const ARint *incx, 
                      const arcomplex<float> *beta, arcomplex<float> *y,
                      const ARint *incy);

  // Double precision complex routines.

  void F77NAME(zdotc)(arcomplex<double> *c, const ARint *n,
                      const arcomplex<double> *cx, const ARint *incx,
                      const arcomplex<double> *cy, const ARint *incy);

  void F77NAME(zdotu)(arcomplex<double> *c, const ARint *n,
                      const arcomplex<double> *cx, const ARint *incx,
                      const arcomplex<double> *cy, const ARint *incy);

  void F77NAME(zaxpy)(const ARint *n, const arcomplex<double> *da,
                      const arcomplex<double> *dx, const ARint *incx,
                      arcomplex<double> *dy, const ARint *incy);

  void F77NAME(zcopy)(const ARint *n, const arcomplex<double> *dx,
                      const ARint *incx, arcomplex<double> *dy,
                      const ARint *incy);

  double  F77NAME(dzasum)(const ARint *n, const arcomplex<double> *dx,
                          const ARint *incx);

  double  F77NAME(dznrm2)(const ARint *n, const arcomplex<double> *dx,
                          const ARint *incx);

  void F77NAME(zdscal)(const ARint *n, const double *da, arcomplex<double> *dx,
                       const ARint *incx);

  void F77NAME(zscal)(const ARint *n, const arcomplex<double> *da,
                      arcomplex<double> *dx, const ARint *incx);

  ARint F77NAME(izamax)(const ARint *n, const arcomplex<double> *dx,
                          const ARint *incx);

  void F77NAME(zswap)(const ARint *n, arcomplex<double> *dx,
                      const ARint *incx, arcomplex<double> *dy, 
                      const ARint *incy);

  void F77NAME(zgemv)(const char* trans, const ARint *m, 
                      const ARint *n, const arcomplex<double> *alpha,
                      const arcomplex<double> *a, const ARint *lda, 
                      const arcomplex<double> *x, const ARint *incx, 
                      const arcomplex<double> *beta, arcomplex<double> *y,
                      const ARint *incy);

  void F77NAME(zgbmv)(const char* trans, const ARint *m, 
                      const ARint *n, const ARint *kl, 
                      const ARint *ku, const arcomplex<double> *alpha,
                      const arcomplex<double> *a, const ARint *lda, 
                      const arcomplex<double> *x, const ARint *incx, 
                      const arcomplex<double> *beta, arcomplex<double> *y,
                      const ARint *incy);

#endif // ARCOMP_H

}
#endif // BLAS1F_H

