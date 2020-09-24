/*
  ARPACK++ v1.2 2/20/2000
  c++ interface to ARPACK code.

  MODULE arpackf.h
  ARPACK FORTRAN routines.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef ARPACKF_H
#define ARPACKF_H

#include "arch.h"

extern "C"
{

// debug "common" statement.

  extern struct { 
    ARint logfil, ndigit, mgetv0;
    ARint msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    ARint mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
    ARint mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
  } F77NAME(debug);


// double precision symmetric routines.

  void F77NAME(dsaupd)(ARint *ido, char *bmat, ARint *n, const char *which,
                       ARint *nev, double *tol, double *resid,
                       ARint *ncv, double *V, ARint *ldv,
                       ARint *iparam, ARint *ipntr, double *workd,
                       double *workl, ARint *lworkl, ARint *info);

  void F77NAME(dseupd)(ARlogical *rvec, char *HowMny, ARlogical *select,
                       double *d, double *Z, ARint *ldz,
                       double *sigma, char *bmat, ARint *n,
                       const char *which, ARint *nev, double *tol,
                       double *resid, ARint *ncv, double *V,
                       ARint *ldv, ARint *iparam, ARint *ipntr,
                       double *workd, double *workl,
                       ARint *lworkl, ARint *info);

// double precision nonsymmetric routines.

  void F77NAME(dnaupd)(ARint *ido, char *bmat, ARint *n, const char *which,
                       ARint *nev, double *tol, double *resid,
                       ARint *ncv, double *V, ARint *ldv,
                       ARint *iparam, ARint *ipntr, double *workd,
                       double *workl, ARint *lworkl, ARint *info);

  void F77NAME(dneupd)(ARlogical *rvec, char *HowMny, ARlogical *select,
                       double *dr, double *di, double *Z,
                       ARint *ldz, double *sigmar,
                       double *sigmai, double *workev,
                       char *bmat, ARint *n, const char *which,
                       ARint *nev, double *tol, double *resid,
                       ARint *ncv, double *V, ARint *ldv,
                       ARint *iparam, ARint *ipntr,
                       double *workd, double *workl,
                       ARint *lworkl, ARint *info);

// single precision symmetric routines.

  void F77NAME(ssaupd)(ARint *ido, char *bmat, ARint *n, const char *which,
                       ARint *nev, float *tol, float *resid,
                       ARint *ncv, float *V, ARint *ldv,
                       ARint *iparam, ARint *ipntr, float *workd,
                       float *workl, ARint *lworkl, ARint *info);

  void F77NAME(sseupd)(ARlogical *rvec, char *HowMny, ARlogical *select,
                       float *d, float *Z, ARint *ldz,
                       float *sigma, char *bmat, ARint *n,
                       const char *which, ARint *nev, float *tol,
                       float *resid, ARint *ncv, float *V,
                       ARint *ldv, ARint *iparam, ARint *ipntr,
                       float *workd, float *workl,
                       ARint *lworkl, ARint *info);

// single precision nonsymmetric routines.

  void F77NAME(snaupd)(ARint *ido, char *bmat, ARint *n, const char *which,
                       ARint *nev, float *tol, float *resid,
                       ARint *ncv, float *V, ARint *ldv,
                       ARint *iparam, ARint *ipntr, float *workd,
                       float *workl, ARint *lworkl, ARint *info);

  void F77NAME(sneupd)(ARlogical *rvec, char *HowMny, ARlogical *select,
                       float *dr, float *di, float *Z,
                       ARint *ldz, float *sigmar,
                       float *sigmai, float *workev, char *bmat,
                       ARint *n, const char *which, ARint *nev,
                       float *tol, float *resid, ARint *ncv,
                       float *V, ARint *ldv, ARint *iparam,
                       ARint *ipntr, float *workd, float *workl,
                       ARint *lworkl, ARint *info);

#ifdef ARCOMP_H

// single precision complex routines.

  void F77NAME(cnaupd)(ARint *ido, char *bmat, ARint *n, const char *which,
                       ARint *nev, float *tol, arcomplex<float> *resid,
                       ARint *ncv, arcomplex<float> *V, ARint *ldv,
                       ARint *iparam, ARint *ipntr, arcomplex<float> *workd,
                       arcomplex<float> *workl, ARint *lworkl,
                       float *rwork, ARint *info);

  void F77NAME(cneupd)(ARlogical *rvec, char *HowMny, ARlogical *select,
                       arcomplex<float> *d, arcomplex<float> *Z, ARint *ldz,
                       arcomplex<float> *sigma, arcomplex<float> *workev,
                       char *bmat, ARint *n, const char *which, ARint *nev,
                       float *tol, arcomplex<float> *resid, ARint *ncv,
                       arcomplex<float> *V, ARint *ldv, ARint *iparam,
                       ARint *ipntr, arcomplex<float> *workd,
                       arcomplex<float> *workl, ARint *lworkl,
                       float *rwork, ARint *info);

// double precision complex routines.

  void F77NAME(znaupd)(ARint *ido, char *bmat, ARint *n, const char *which,
                       ARint *nev, double *tol, arcomplex<double> *resid,
                       ARint *ncv, arcomplex<double> *V, ARint *ldv,
                       ARint *iparam, ARint *ipntr, arcomplex<double> *workd,
                       arcomplex<double> *workl, ARint *lworkl,
                       double *rwork, ARint *info);

  void F77NAME(zneupd)(ARlogical *rvec, char *HowMny, ARlogical *select,
                       arcomplex<double> *d, arcomplex<double> *Z, ARint *ldz,
                       arcomplex<double> *sigma, arcomplex<double> *workev,
                       char *bmat, ARint *n, const char *which, ARint *nev,
                       double *tol, arcomplex<double> *resid, ARint *ncv,
                       arcomplex<double> *V, ARint *ldv, ARint *iparam,
                       ARint *ipntr, arcomplex<double> *workd,
                       arcomplex<double> *workl, ARint *lworkl,
                       double *rwork, ARint *info);

}

#endif // ARCOMP_H

#endif // ARPACKF_H
