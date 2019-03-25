#ifndef HEADER_myblas
#define HEADER_myblas

/* ************************************************************************ */
/* BLAS function interface with local and external loadable versions        */
/* Author:  Kjell Eikland                                                   */
/* Version: Initial version spring 2004                                     */
/* Licence: LGPL                                                            */
/* ************************************************************************ */
/* Changes: 19 September 2004   Moved function pointer variable             */
/*                              declarations from myblas.h to myblas.c      */
/*                              to avoid linker problems with the Mac.      */
/*          20 April 2005       Modified all double types to REAL to self-  */
/*                              adjust to global settings.  Note that BLAS  */
/*                              as of now does not have double double.      */
/*          11/11/2014          Modified for R package.                     */
/* ************************************************************************ */

#define BLAS_BASE         1
#include "commonlib.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ************************************************************************ */
/* BLAS functions                                                           */
/* ************************************************************************ */

void init_BLAS(void);
MYBOOL is_nativeBLAS(void);
MYBOOL load_BLAS(char *libname);
MYBOOL unload_BLAS(void);

/* ************************************************************************ */
/* User-callable BLAS definitions (C base 1)                                */
/* ************************************************************************ */

void lps_dscal ( int n, LPSREAL da,  LPSREAL *dx, int incx );
void lps_dcopy ( int n, LPSREAL *dx, int incx, LPSREAL *dy, int incy );
void lps_daxpy ( int n, LPSREAL da,  LPSREAL *dx, int incx,   LPSREAL *dy, int incy );
int  lps_idamax ( int n, LPSREAL *x,  int is );


#ifdef __cplusplus
}
#endif

#endif
