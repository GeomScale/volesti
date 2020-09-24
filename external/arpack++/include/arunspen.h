/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUNSPen.h.
   Arpack++ class ARumNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUNSPEN_H
#define ARUNSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "umfpackc.h"
#include "arunsmat.h"


template<class ARTYPE, class ARFLOAT>
class ARumNonSymPencil
{

 protected:

  char                               part;
  ARumNonSymMatrix<ARTYPE, ARFLOAT>* A;
  ARumNonSymMatrix<ARTYPE, ARFLOAT>* B;
  ARumNonSymMatrix<ARTYPE, ARFLOAT>  AsB;
#ifdef ARCOMP_H
  ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> AsBc;
#endif

  virtual void Copy(const ARumNonSymPencil& other);

  void SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
                   int yind[], int ny, ARTYPE z[], int zind[], int& nz);

#ifdef ARCOMP_H
  void SparseSaxpy(arcomplex<ARFLOAT> a, ARFLOAT x[], int xind[], int nx,
                   ARFLOAT y[], int yind[], int ny,
                   arcomplex<ARFLOAT> z[], int zind[], int& nz);
#endif

  void SubtractAsB(ARTYPE sigma);

#ifdef ARCOMP_H
  void SubtractAsB(ARFLOAT sigmaR, ARFLOAT sigmaI);
#endif

 public:

#ifdef ARCOMP_H
  bool IsFactored() { return (AsB.IsFactored()||AsBc.IsFactored()); }
#else
  bool IsFactored() { return AsB.IsFactored(); }
#endif

  bool IsSymmetric() { return AsB.IsSymmetric(); }

  void FactorAsB(ARTYPE sigma);

#ifdef ARCOMP_H
  void FactorAsB(ARFLOAT sigmaR, ARFLOAT sigmaI, char partp = 'R');
#endif

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

#ifdef ARCOMP_H
  void MultInvAsBv(arcomplex<ARFLOAT>* v, arcomplex<ARFLOAT>* w);
#endif

  void MultInvAsBv(ARFLOAT* v, ARFLOAT* w);

  void DefineMatrices(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                      ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp);

  ARumNonSymPencil() { part = 'N'; }
  // Short constructor that does nothing.

  ARumNonSymPencil(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                   ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp);
  // Long constructor.

  ARumNonSymPencil(const ARumNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumNonSymPencil() { }
  // Destructor.

  ARumNonSymPencil& operator=(const ARumNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymPencil<ARTYPE, ARFLOAT>::
Copy(const ARumNonSymPencil<ARTYPE, ARFLOAT>& other)
{

  part     = other.part;
  A        = other.A;
  B        = other.B;
  AsB      = other.AsB;
#ifdef ARCOMP_H
  AsBc     = other.AsBc;
#endif

} // Copy.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
            int yind[], int ny, ARTYPE z[], int zind[], int& nz)
// A strongly sequential (and inefficient) sparse saxpy algorithm.
{

  int ix, iy;

  nz = 0;
  if ((nx == 0) || (a == (ARTYPE)0)) {
    copy(ny,y,1,z,1);
    for (iy=0; iy!=ny; iy++) zind[iy] = yind[iy];
    nz = ny;
    return;
  }
  if (ny == 0) {
    copy(nx,x,1,z,1);
    scal(nx,a,z,1);
    for (ix=0; ix!=nx; ix++) zind[ix] = xind[ix];
    nz = nx;
    return;
  }
  ix = 0;
  iy = 0;
  while (true) {
    if (xind[ix] == yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++]+y[iy++];
      if ((ix == nx)||(iy == ny)) break;
    }
    else if (xind[ix] < yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++];
      if (ix == nx) break;
    }
    else {
      zind[nz] = yind[iy];
      z[nz++]  = y[iy++];
      if (iy == ny) break;
    }
  }
  while (iy < ny) {
    zind[nz] = yind[iy];
    z[nz++]  = y[iy++];
  }
  while (ix < nx) {
    zind[nz] = xind[ix];
    z[nz++]  = x[ix++];
  }

} // SparseSaxpy (ARTYPE).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
SparseSaxpy(arcomplex<ARFLOAT> a, ARFLOAT x[], int xind[], int nx,
            ARFLOAT y[], int yind[], int ny, 
            arcomplex<ARFLOAT> z[], int zind[], int& nz)
// A strongly sequential (and inefficient) sparse saxpy algorithm.
{

  int ix, iy;

  nz = 0;
  if ((nx == 0) || (a == arcomplex<ARFLOAT>(0.0,0.0))) {
    for (iy=0; iy!=ny; iy++) {
      z[iy]    = arcomplex<ARFLOAT>(y[iy],0.0);
      zind[iy] = yind[iy];
    }
    nz = ny;
    return;
  }
  if (ny == 0) {
    for (ix=0; ix!=ny; ix++) {
      z[ix]    = a*arcomplex<ARFLOAT>(x[ix],0.0);
      zind[ix] = xind[ix];
    }
    nz = nx;
    return;
  }
  ix = 0;
  iy = 0;
  while (true) {
    if (xind[ix] == yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++]+y[iy++];
      if ((ix == nx)||(iy == ny)) break;
    }
    else if (xind[ix] < yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++];
      if (ix == nx) break;
    }
    else {
      zind[nz] = yind[iy];
      z[nz++]  = arcomplex<ARFLOAT>(y[iy++], 0.0);
      if (iy == ny) break;
    }
  }
  while (iy < ny) {
    zind[nz] = yind[iy];
    z[nz++]  = arcomplex<ARFLOAT>(y[iy++], 0.0);
  }
  while (ix < nx) {
    zind[nz] = xind[ix];
    z[nz++]  = arcomplex<ARFLOAT>(x[ix++], 0.0);
  }

} // SparseSaxpy (arcomplex<ARFLOAT>).
#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::SubtractAsB(ARTYPE sigma)
{

  int i, acol, bcol, asbcol, scol;

  // Subtracting sigma*B from A.

  AsB.index[0] = 0;
  asbcol       = 0;

  for (i=0; i!=AsB.n; i++) {
    bcol = B->pcol[i];
    acol = A->pcol[i];
    SparseSaxpy(-sigma, &B->a[bcol], &B->irow[bcol], B->pcol[i+1]-bcol,
                &A->a[acol], &A->irow[acol], A->pcol[i+1]-acol,
                &AsB.value[asbcol], &AsB.index[asbcol+AsB.n+1], scol);
    asbcol += scol;
    AsB.index[i+1] = asbcol;
  }

  AsB.nnz = asbcol; 

  // Adding one to all elements of vector index
  // because the decomposition function was written in FORTRAN.

  for (i=0; i<=AsB.n+AsB.nnz; i++) AsB.index[i]++;

} // SubtractAsB (ARTYPE shift).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
SubtractAsB(ARFLOAT sigmaR, ARFLOAT sigmaI)
{

  int                i, acol, bcol, asbcol, scol;
  arcomplex<ARFLOAT> sigma;

  // Subtracting sigma*B from A.

  sigma         = arcomplex<ARFLOAT>(sigmaR, sigmaI);
  AsBc.index[0] = 0;
  asbcol        = 0;

  for (i=0; i!=AsBc.n; i++) {
    bcol = B->pcol[i];
    acol = A->pcol[i];
    SparseSaxpy(-sigma, &B->a[bcol], &B->irow[bcol], B->pcol[i+1]-bcol,
                &A->a[acol], &A->irow[acol], A->pcol[i+1]-acol,
                &AsBc.value[asbcol], &AsBc.index[asbcol+AsBc.n+1], scol);
    asbcol += scol;
    AsBc.index[i+1] = asbcol;
  }

  AsBc.nnz = asbcol;

  // Adding one to all elements of vector index
  // because the decomposition function was written in FORTRAN.

  for (i=0; i<=AsBc.n+AsBc.nnz; i++) AsBc.index[i]++;

} // SubtractAsB (arcomplex<ARFLOAT> shift).
#endif // ARCOMP_H


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {

    int fillin = A->fillin > B->fillin ? A->fillin : B->fillin;
    AsB.DefineMatrix(A->ncols(), A->nzeros(), A->a, A->irow, 
                     A->pcol, A->threshold, fillin, 
                     (A->IsSymmetric() && B->IsSymmetric()), 
                     A->icntl[3], false);
    AsB.nnz = A->nzeros()+B->nzeros(); // temporary value.

  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure(); // AsB.nnz must be set to A->nzeros()+B->nzeros().

  // Subtracting sigma*B from A and storing the result on AsB.

  SubtractAsB(sigma);

  // Decomposing AsB.

  um2fa(AsB.n, AsB.index[AsB.n], 0, false, AsB.lvalue, AsB.lindex, AsB.value,
        AsB.index, AsB.keep, AsB.cntl, AsB.icntl, AsB.info, AsB.rinfo);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (ARTYPE shift).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
FactorAsB(ARFLOAT sigmaR, ARFLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsBc.IsDefined()) {

    part        = partp;
    int  fillin = A->fillin > B->fillin ? A->fillin : B->fillin;
    AsBc.DefineMatrix(A->ncols(), A->nzeros(), 0, 0,
                      A->pcol, A->threshold, fillin, 
                      (A->IsSymmetric() && B->IsSymmetric()), 
                      A->icntl[3], false);
    AsBc.nnz    = A->nzeros()+B->nzeros(); // temporary value.

  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsBc.CreateStructure(); // AsBc.nnz must be set to A->nzeros()+B->nzeros().

  // Subtracting sigma*B from A and storing the result on AsBc.

  SubtractAsB(sigmaR, sigmaI);

  // Decomposing AsB.

  um2fa(AsBc.n, AsBc.index[AsBc.n], 0, false, AsBc.lvalue, AsBc.lindex,
        AsBc.value, AsBc.index, AsBc.keep, AsBc.cntl, AsBc.icntl, 
        AsBc.info, AsBc.rinfo);

  // Handling errors.

  AsBc.ThrowError();

  AsBc.factored = true;

} // FactorAsB (arcomplex<ARFLOAT> shift).
#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H

template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
MultInvAsBv(arcomplex<ARFLOAT>* v, arcomplex<ARFLOAT>* w)
{

  AsB.MultInvv((ARTYPE*)v,(ARTYPE*)w);

} // MultInvAsBv (arcomplex<ARFLOAT>).

#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::MultInvAsBv(ARFLOAT* v, ARFLOAT* w)
{

  if (part == 'N') {    // shift is real.

    AsB.MultInvv((ARTYPE*)v,(ARTYPE*)w);

  }
  else {                // shift is complex.

#ifdef ARCOMP_H

    int                i;
    arcomplex<ARFLOAT> *tv, *tw;

    tv = new arcomplex<ARFLOAT>[AsBc.ncols()];
    tw = new arcomplex<ARFLOAT>[AsBc.ncols()];

    for (i=0; i!=AsBc.ncols(); i++) tv[i] = arcomplex<ARFLOAT>(v[i], 0.0);

    AsBc.MultInvv(tv, tw);

    if (part=='I') {
      for (i=0; i!=AsBc.ncols(); i++) w[i] = imag(tw[i]);
    }
    else {
      for (i=0; i!=AsBc.ncols(); i++) w[i] = real(tw[i]);
    }

    delete[] tv;
    delete[] tw;

#endif // ARCOMP_H.

  }

} // MultInvAsBv (ARFLOAT).


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymPencil<ARTYPE, ARFLOAT>::
DefineMatrices(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
               ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  A = &Ap;
  B = &Bp;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARumNonSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE, class ARFLOAT>
inline ARumNonSymPencil<ARTYPE, ARFLOAT>::
ARumNonSymPencil(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                 ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE, class ARFLOAT>
ARumNonSymPencil<ARTYPE, ARFLOAT>& ARumNonSymPencil<ARTYPE, ARFLOAT>::
operator=(const ARumNonSymPencil<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUNSPEN_H
