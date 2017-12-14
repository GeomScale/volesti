/*
FILE DETAILS
description: header file needed for LUmod functions 
project: MPT 3.0
filename: lumod_dense.h
author:  Michael Saunders, Systems Optimization Laboratory, Department of EESOR, Stanford
University. 

REVISION HISTORY:
date: Nov 29, 2010
revised by: Martin Herceg, Automatic Control Laboratory, ETH Zurich, 2010
details: Added declarations of TINYNUMBER, MACHINEPREC and replaced function arguments int with
ptrdiff_t to be compatible with Matlab. Added declarations of customized Blas functions. Also added a license text.
 
LICENSE:

     BSD License for SOL Numerical Software
    http://www.opensource.org/licenses/bsd-license.php

This notice applies to the software packages

    cgLanczos, cgls, lsmr, lsqr, lumod, lusol, minres, minres-qlp, pdco, symmlq

made available at
    http://www.stanford.edu/group/SOL/software.html
by the
    Systems Optimization Laboratory (SOL)
    Dept of Management Science and Engineering
    Stanford University,
    Stanford, CA 94305-4026, USA


Copyright (c) 2010, Systems Optimization Laboratory
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of Stanford University nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef LUMOD_DENSE__H
#define LUMOD_DENSE__H

#define TINYNUMBER     1.0e-04
#define MACHINEPREC   2.22e-16
#define subvec(item) (item -1)

/* LU functions */
void LUmod( int mode, ptrdiff_t maxmod, ptrdiff_t n, ptrdiff_t krow, ptrdiff_t kcol,
             double *L, double *U, double *y, double *z, double *w );
void Lprod ( int mode, ptrdiff_t maxmod, ptrdiff_t n,
             double *L, double *y, double *z );
void LUforw ( ptrdiff_t first, ptrdiff_t last, ptrdiff_t n, ptrdiff_t nu, ptrdiff_t maxmod,
              double eps, double *L, double *U, double *y );
void LUback ( ptrdiff_t first, ptrdiff_t *last, ptrdiff_t n, ptrdiff_t nu,
              ptrdiff_t maxmod, double eps,
              double *L, double *U, double *y, double *z );
void Usolve ( int mode, ptrdiff_t maxmod, ptrdiff_t n,
              double *U, double *y );
void Lprod ( int mode, ptrdiff_t maxmod, ptrdiff_t n,
             double *L, double *y, double *z );
void elm ( ptrdiff_t first, ptrdiff_t last, double *x, double *y, double cs, double sn );
void elmgen ( double *x, double *y, double eps, double *cs, double *sn );

/* customized BLAS functions */
void LUdaxpy ( ptrdiff_t n, double da, double *dx, ptrdiff_t incx,double *dy, ptrdiff_t incy);
void LUdcopy ( ptrdiff_t n, double *dx, ptrdiff_t incx, double *dy, ptrdiff_t incy);
double LUddot ( ptrdiff_t n, double *dx, ptrdiff_t incx, double *dy, ptrdiff_t incy);
void LUdswap ( ptrdiff_t n, double *dx, ptrdiff_t incx, double *dy, ptrdiff_t incy );

#endif
