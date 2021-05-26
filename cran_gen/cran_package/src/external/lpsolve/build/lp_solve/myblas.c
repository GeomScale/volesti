
#include <stdlib.h>
#include <stdio.h>
/*#include <memory.h>*/
#include <string.h>
#include <math.h>
#include "myblas.h"
#include "R_ext/BLAS.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif

/* ************************************************************************ */
/* Initialize BLAS interfacing routines                                     */
/* ************************************************************************ */

MYBOOL mustinitBLAS = TRUE;
#ifdef WIN32
  HINSTANCE hBLAS = NULL;
#else
  void      *hBLAS = NULL;
#endif


/* ************************************************************************ */
/* Define the BLAS interfacing routines                                     */
/* ************************************************************************ */

void init_BLAS(void)
{
  if(mustinitBLAS) {
    load_BLAS(NULL);
    mustinitBLAS = FALSE;
  }
}

MYBOOL is_nativeBLAS(void)
{
#ifdef LoadableBlasLib
  return( (MYBOOL) (hBLAS == NULL) );
#else
  return( TRUE );
#endif
}

MYBOOL load_BLAS(char *libname)
{
  MYBOOL result = TRUE;
  return( result );
}

MYBOOL unload_BLAS(void)
{
  return( load_BLAS(NULL) );
}


/* ************************************************************************ */
/* Now define the unoptimized local BLAS functions                          */
/* ************************************************************************ */

void lps_daxpy( int n, LPSREAL da, LPSREAL *dx, int incx, LPSREAL *dy, int incy)
{
  dx++;
  dy++;
  F77_CALL(daxpy) ( &n, &da, dx, &incx, dy, &incy);
}


/* ************************************************************************ */

void lps_dcopy( int n, LPSREAL *dx, int incx, LPSREAL *dy, int incy)
{
  dx++;
  dy++;
  F77_CALL(dcopy) ( &n, dx, &incx, dy, &incy);
}


/* ************************************************************************ */

void lps_dscal (int n, LPSREAL da, LPSREAL *dx, int incx)
{
  dx++;
  F77_CALL(dscal) (&n, &da, dx, &incx);
}


/* ************************************************************************ */

int lps_idamax( int n, LPSREAL *x, int is )
{
  x++;
  return ( F77_CALL(idamax)( &n, x, &is ) );
}





