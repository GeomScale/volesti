/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE SuperLUc.h.
   Interface to SuperLU routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SUPERLUC_H
#define SUPERLUC_H

#include "arch.h"
#include "arlspdef.h"
#include "arlsupm.h"
#include "arlcomp.h"

// gstrf.

inline void gstrf(superlu_options_t *options, SuperMatrix *A,
        int relax, int panel_size, int *etree, void *work, int lwork,
        int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
        SuperLUStat_t *stat, int *info)
{
  if (A->Dtype == SLU_D) {       // calling the double precision routine.
    dGlobalLU_t Glu;
    dgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
  }
  else if (A->Dtype == SLU_S) {  // calling the single precision routine.
    sGlobalLU_t Glu;
    sgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
  }
  else if (A->Dtype == SLU_Z) {  // calling the double precision complex routine.
#ifdef ARCOMP_H
    zGlobalLU_t Glu;
    zgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
#endif
  }
  else {                      // calling the single precision complex routine.
#ifdef ARCOMP_H
    cGlobalLU_t Glu;
    cgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
#endif
  }

} // gstrf.


inline void gstrs(trans_t trans, SuperMatrix *L, SuperMatrix *U,
	          int *perm_c, int *perm_r, SuperMatrix *B, SuperLUStat_t* stat, int *info)
{

  if (L->Dtype == SLU_D) {       // calling the double precision routine.
    dgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
  }
  else if (L->Dtype == SLU_S) {  // calling the single precision routine.
    sgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
  }
  else if (L->Dtype == SLU_Z) {  // calling the double precision complex routine.
#ifdef ARCOMP_H
    zgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
#endif
  }
  else {                      // calling the single precision complex routine.
#ifdef ARCOMP_H
    cgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
#endif
  }

} // gstrs.


// Create_CompCol_Matrix.

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  double* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  dCreate_CompCol_Matrix(A,m,n,nnz,a,irow,pcol,S,SLU_D,M);

} // Create_CompCol_Matrix (double).

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  float* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  sCreate_CompCol_Matrix(A,m,n,nnz,a,irow,pcol,S,SLU_S,M);

} // Create_CompCol_Matrix (float).

#ifdef ARCOMP_H

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  arcomplex<double>* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  zCreate_CompCol_Matrix(A,m,n,nnz,(ldcomplex*)a,irow,pcol,S,SLU_Z,M);

} // Create_CompCol_Matrix (complex<double>).

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  arcomplex<float>* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  cCreate_CompCol_Matrix(A,m,n,nnz,(lscomplex*)a,irow,pcol,S,SLU_C,M);

} // Create_CompCol_Matrix (complex<float>).

#endif // ARCOMP_H.


// Create_Dense_Matrix.

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, double* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  dCreate_Dense_Matrix(A,m,n,x,ldx,S,SLU_D,M);

} // Create_Dense_Matrix (double).

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, float* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  sCreate_Dense_Matrix(A,m,n,x,ldx,S,SLU_S,M);

} // Create_Dense_Matrix (float).

#ifdef ARCOMP_H

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, arcomplex<double>* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  zCreate_Dense_Matrix(A,m,n,(ldcomplex*)x,ldx,S,SLU_Z,M);

} // Create_Dense_Matrix (complex<double>).

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, arcomplex<float>* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  cCreate_Dense_Matrix(A,m,n,(lscomplex*)x,ldx,S,SLU_C,M);

} // Create_Dense_Matrix (complex<float>).

#endif // ARCOMP_H.

#endif // SUPERLUC_H
