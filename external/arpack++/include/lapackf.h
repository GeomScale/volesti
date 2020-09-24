/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE lapackf.h.
   Interface to LAPACK FORTRAN routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LAPACKF_H
#define LAPACKF_H

#include "arch.h"

extern "C"
{

  // Single precision real routines.

  float F77NAME(slapy2)(const float *x, const float *y);

  void F77NAME(slacpy)(const char* uplo, const ARint *m, const ARint *n,
                       const float *a, const ARint *lda, float *b, 
                       const ARint *ldb);

  void F77NAME(sgttrf)(const ARint *n, float *dl, float *d, float *du,
                       float *du2, ARint *ipiv, ARint *info);

  void F77NAME(sgbtrf)(const ARint *m, const ARint *n, const ARint *kl,
                       const ARint *ku, float *ab, const ARint *ldab,
                       ARint *ipiv, ARint *info);

  void F77NAME(sgetrf)(const ARint *m, const ARint *n, float *A,
                       const ARint *lda, ARint *ipiv, ARint *info);

  void F77NAME(sgttrs)(const char* trans, const ARint *n,
                       const ARint *nrhs, const float *dl,
                       const float *d, const float *du,
                       const float *du2, const ARint *ipiv,
                       float* b, const ARint *ldb, ARint *info);

  void F77NAME(sgbtrs)(const char* trans, const ARint *n, 
                       const ARint *kl, const ARint *ku, 
                       const ARint *nrhs, const float *ab, 
                       const ARint *ldab, const ARint *ipiv, 
                       float *b, const ARint *ldb, ARint *info);

  void F77NAME(sgetrs)(const char* trans, const ARint *n,
                       const ARint *nrhs, const float *A,
                       const ARint *lda, const ARint *ipiv,
                       float* b, const ARint *ldb, ARint *info);

  void F77NAME(spttrf)(const ARint *n, float *d, float *e, ARint *info);

  void F77NAME(spttrs)(const ARint *n, const ARint *nrhs,
                       const float *d, const float *e, float *b,
                       const ARint *ldb, ARint *info);

  void F77NAME(ssptrf)(const char* trans, const ARint *n, 
                       float *ap, ARint *ipiv, ARint *info);

  void F77NAME(ssptrs)(const char* trans, const ARint *n, 
                       const ARint *nrhs, float *ap, ARint *ipiv, 
                       float *b, const ARint *ldb, ARint *info);

  // Double precision real routines.

  double F77NAME(dlapy2)(const double *x, const double *y);

  void F77NAME(dlacpy)(const char* uplo, const ARint *m, const ARint *n,
                       const double *a, const ARint *lda, double *b, 
                       const ARint *ldb);

  void F77NAME(dgttrf)(const ARint *n, double *dl, double *d, double *du,
                       double *du2, ARint *ipiv, ARint *info);

  void F77NAME(dgbtrf)(const ARint *m, const ARint *n, const ARint *kl,
                       const ARint *ku, double *ab, const ARint *ldab,
                       ARint *ipiv, ARint *info);

  void F77NAME(dgetrf)(const ARint *m, const ARint *n, double *A,
                       const ARint *lda, ARint *ipiv, ARint *info);

  void F77NAME(dgttrs)(const char* trans, const ARint *n,
                       const ARint *nrhs, const double *dl,
                       const double *d, const double *du,
                       const double *du2, const ARint *ipiv,
                       double* b, const ARint *ldb, ARint *info);

  void F77NAME(dgbtrs)(const char* trans, const ARint *n, 
                       const ARint *kl, const ARint *ku, 
                       const ARint *nrhs, const double *ab, 
                       const ARint *ldab, const ARint *ipiv, 
                       double *b, const ARint *ldb, ARint *info);

  void F77NAME(dgetrs)(const char* trans, const ARint *n,
                       const ARint *nrhs, const double *A,
                       const ARint *lda, const ARint *ipiv,
                       double* b, const ARint *ldb, ARint *info);

  void F77NAME(dpttrf)(const ARint *n, double *d, double *e, ARint *info);

  void F77NAME(dpttrs)(const ARint *n, const ARint *nrhs,
                       const double *d, const double *e, double *b,
                       const ARint *ldb, ARint *info);

  void F77NAME(dsptrf)(const char* trans, const ARint *n, 
                       double *ap, ARint *ipiv, ARint *info);

  void F77NAME(dsptrs)(const char* trans, const ARint *n, 
                       const ARint *nrhs, double *ap, ARint *ipiv, 
                       double *b, const ARint *ldb, ARint *info);

#ifdef ARCOMP_H

  // Single precision complex routines.

  void F77NAME(clacpy)(const char* uplo, const ARint *m, const ARint *n,
                       const arcomplex<float> *a, const ARint *lda, 
                       arcomplex<float> *b, const ARint *ldb);

  void F77NAME(cgttrf)(const ARint *n, arcomplex<float> *dl,
                       arcomplex<float> *d, arcomplex<float> *du,
                       arcomplex<float> *du2, ARint *ipiv,
                       ARint *info);

  void F77NAME(cgbtrf)(const ARint *m, const ARint *n, const ARint *kl,
                       const ARint *ku, arcomplex<float> *ab, 
                       const ARint *ldab, ARint *ipiv, ARint *info);

  void F77NAME(cgetrf)(const ARint *m, const ARint *n, arcomplex<float> *A,
                       const ARint *lda, ARint *ipiv, ARint *info);

  void F77NAME(cgttrs)(const char *trans, const ARint *n,
                       const ARint *nrhs, const arcomplex<float> *dl,
                       const arcomplex<float> *d, const arcomplex<float> *du,
                       const arcomplex<float> *du2, const ARint *ipiv,
                       arcomplex<float>* b, const ARint *ldb,
                       ARint *info);

  void F77NAME(cgbtrs)(const char* trans, const ARint *n, 
                       const ARint *kl, const ARint *ku, 
                       const ARint *nrhs, const arcomplex<float> *ab, 
                       const ARint *ldab, const ARint *ipiv, 
                       arcomplex<float> *b, const ARint *ldb, 
                       ARint *info);

  void F77NAME(cgetrs)(const char* trans, const ARint *n,
                       const ARint *nrhs, const arcomplex<float> *A,
                       const ARint *lda, const ARint *ipiv,
                       arcomplex<float>* b, const ARint *ldb, ARint *info);

  // Double precision complex routines.

  void F77NAME(zlacpy)(const char* uplo, const ARint *m, const ARint *n,
                       const arcomplex<double> *a, const ARint *lda, 
                       arcomplex<double> *b, const ARint *ldb);

  void F77NAME(zgttrf)(const ARint *n, arcomplex<double> *dl,
                       arcomplex<double> *d, arcomplex<double> *du,
                       arcomplex<double> *du2, ARint *ipiv,
                       ARint *info);

  void F77NAME(zgbtrf)(const ARint *m, const ARint *n, const ARint *kl,
                       const ARint *ku, arcomplex<double> *ab, 
                       const ARint *ldab, ARint *ipiv, ARint *info);

  void F77NAME(zgetrf)(const ARint *m, const ARint *n, arcomplex<double> *A,
                       const ARint *lda, ARint *ipiv, ARint *info);

  void F77NAME(zgttrs)(const char *trans, const ARint *n,
                       const ARint *nrhs, const arcomplex<double> *dl,
                       const arcomplex<double> *d, const arcomplex<double> *du,
                       const arcomplex<double> *du2, const ARint *ipiv,
                       arcomplex<double>* b, const ARint *ldb,
                       ARint *info);

  void F77NAME(zgbtrs)(const char* trans, const ARint *n, 
                       const ARint *kl, const ARint *ku, 
                       const ARint *nrhs, const arcomplex<double> *ab, 
                       const ARint *ldab, const ARint *ipiv, 
                       arcomplex<double> *b, const ARint *ldb, 
                       ARint *info);

  void F77NAME(zgetrs)(const char* trans, const ARint *n,
                       const ARint *nrhs, const arcomplex<double> *A,
                       const ARint *lda, const ARint *ipiv,
                       arcomplex<double>* b, const ARint *ldb, ARint *info);

#endif // ARCOMP_H

  void F77NAME(second)(const float *T);

}
#endif // LAPACKF_H





