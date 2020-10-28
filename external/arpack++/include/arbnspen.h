/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBNSPen.h.
   Arpack++ class ARbdNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBNSPEN_H
#define ARBNSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arbnsmat.h"


template<class ARTYPE, class ARFLOAT>
class ARbdNonSymPencil
{

 protected:

  char                               part;
  ARbdNonSymMatrix<ARTYPE, ARFLOAT>* A;
  ARbdNonSymMatrix<ARTYPE, ARFLOAT>* B;
  ARbdNonSymMatrix<ARTYPE, ARFLOAT>  AsB;
#ifdef ARCOMP_H
  ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> AsBc;
#endif

  int max(int a, int b) { return (a>b)?a:b; }

  int min(int a, int b) { return (a<b)?a:b; }

  void ComplexCopy(int n, ARFLOAT dx[], int incx, 
                   arcomplex<ARFLOAT> dy[], int incy);

  void ComplexAxpy(int n, arcomplex<ARFLOAT> da, ARTYPE dx[], 
                   int incx, arcomplex<ARFLOAT> dy[], int incy);

  virtual void Copy(const ARbdNonSymPencil& other);

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

  void DefineMatrices(ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                      ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Bp);

  ARbdNonSymPencil() { part = 'N'; }
  // Short constructor that does nothing.

  ARbdNonSymPencil(ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                   ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Bp);
  // Long constructor.

  ARbdNonSymPencil(const ARbdNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdNonSymPencil() { }
  // Destructor.

  ARbdNonSymPencil& operator=(const ARbdNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARbdNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
inline void ARbdNonSymPencil<ARTYPE, ARFLOAT>::
Copy(const ARbdNonSymPencil<ARTYPE, ARFLOAT>& other)
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
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::
ComplexCopy(int n, ARFLOAT dx[], int incx, arcomplex<ARFLOAT> dy[], int incy)
{

  for (int ix=0, iy=0; ix<(n*incx); ix+=incx, iy+=incy) {
    dy[iy] = arcomplex<ARFLOAT>(dx[ix], 0.0);
  }

} // ComplexCopy.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::
ComplexAxpy(int n, arcomplex<ARFLOAT> da, ARTYPE dx[], int incx,
            arcomplex<ARFLOAT> dy[], int incy)
{

  for (int ix=0, iy=0; ix<(n*incx); ix+=incx, iy+=incy) {
    dy[iy] += da*dx[ix];
  }

} // ComplexAxpy.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::SubtractAsB(ARTYPE sigma)
{

  int    i, inca, incb, minL, minU, begB, begAsB;
  ARTYPE negsig;

  inca   = A->ndiagL+A->ndiagU+1;
  incb   = B->ndiagL+B->ndiagU+1;
  negsig = -sigma;

  // Expanding A.

  begAsB = AsB.ndiagL+AsB.ndiagU-A->ndiagU;
  for (i = 0; i < inca; i++) {
    copy(AsB.n, &A->A[i], inca, &AsB.Ainv[begAsB+i], AsB.lda);
  }

  // Transferring part of B (*(-sigma)) if AsB.ndiagU > A->ndiagU.

  if (A->ndiagU < AsB.ndiagU) {
    for (i = 0; i < AsB.ndiagU-A->ndiagU; i++) {
      copy(AsB.n, &B->A[i], incb, &AsB.Ainv[AsB.ndiagL+i], AsB.lda);
      scal(AsB.n, negsig, &AsB.Ainv[AsB.ndiagL+i], AsB.lda);
    }
  } 

  // Subtracting sigma*B from A.

  minL   = min(A->ndiagL, B->ndiagL);
  minU   = min(A->ndiagU, B->ndiagU);
  begB   = B->ndiagU-minU;
  begAsB = AsB.ndiagL+AsB.ndiagU-minU;

  for (i = 0; i < minL+minU+1; i++) {
    axpy(AsB.n, -sigma, &B->A[begB+i], incb, &AsB.Ainv[begAsB+i], AsB.lda);  
  }

  // Transferring part of B (*(-sigma)) if AsB.ndiagL > A->ndiagL.

  if (A->ndiagL < AsB.ndiagL) {
    begB   = B->ndiagU+1+minL;
    begAsB = AsB.ndiagL+AsB.ndiagU+1+minL;
    for (i = 0; i < AsB.ndiagL-A->ndiagL; i++) {
      copy(AsB.n, &B->A[begB+i], incb, &AsB.Ainv[begAsB+i], AsB.lda);
      scal(AsB.n, negsig, &AsB.Ainv[begAsB+i], AsB.lda);
    }
  } 

} // SubtractAsB (ARTYPE shift).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::
SubtractAsB(ARFLOAT sigmaR, ARFLOAT sigmaI)
{

  int                i, inca, incb, minL, minU, begB, begAsB;
  arcomplex<ARFLOAT> sigma;

  inca  = A->ndiagL+A->ndiagU+1;
  incb  = B->ndiagL+B->ndiagU+1;
  sigma = arcomplex<ARFLOAT>(sigmaR, sigmaI);

  // Expanding A.

  begAsB = AsBc.ndiagL+AsBc.ndiagU-A->ndiagU;
  for (i = 0; i < inca; i++) {
    ComplexCopy(AsBc.n,(ARFLOAT*)(&A->A[i]),inca,&AsBc.Ainv[begAsB+i],AsBc.lda);
  }

  // Transferring part of B (*(-sigma)) if AsBc.ndiagU > A->ndiagU.

  if (A->ndiagU < AsBc.ndiagU) {
    for (i = 0; i < AsBc.ndiagU-A->ndiagU; i++) {
      ComplexCopy(AsBc.n, (ARFLOAT*)(&B->A[i]), incb,
                  &AsBc.Ainv[AsBc.ndiagL+i], AsBc.lda);
      scal(AsBc.n, -sigma, &AsBc.Ainv[AsBc.ndiagL+i], AsBc.lda);
    }
  } 

  // Subtracting sigma*B from A.

  minL   = min(A->ndiagL, B->ndiagL);
  minU   = min(A->ndiagU, B->ndiagU);
  begB   = B->ndiagU-minU;
  begAsB = AsBc.ndiagL+AsBc.ndiagU-minU;

  for (i = 0; i < minL+minU+1; i++) {
    ComplexAxpy(AsBc.n, -sigma, &B->A[begB+i], incb, 
                &AsBc.Ainv[begAsB+i], AsBc.lda);  
  }

  // Transferring part of B (*(-sigma)) if AsBc.ndiagL > A->ndiagL.

  if (A->ndiagL < AsBc.ndiagL) {
    begB   = B->ndiagU+1+minL;
    begAsB = AsBc.ndiagL+AsBc.ndiagU+1+minL;
    for (i = 0; i < AsBc.ndiagL-A->ndiagL; i++) {
      ComplexCopy(AsBc.n, (ARFLOAT*)(&B->A[begB+i]), incb,
                  &AsBc.Ainv[begAsB+i], AsBc.lda);
      scal(AsBc.n, -sigma, &AsBc.Ainv[begAsB+i], AsBc.lda);
    }
  }

} // SubtractAsB (arcomplex<ARFLOAT> shift).
#endif // ARCOMP_H


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Copying A to AsB if sigma = 0.

  if (sigma == (ARTYPE)0) {

    AsB = *A;
    if (!AsB.IsFactored()) AsB.FactorA();
    return;

  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {
    AsB.DefineMatrix(A->ncols(), max(A->ndiagL, B->ndiagL),
                     max(A->ndiagU, B->ndiagU), A->A);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure(); 

  // Subtracting sigma*B from A and storing the result on AsB.

  SubtractAsB(sigma);

  // Decomposing AsB.

  gbtrf(AsB.n, AsB.n, AsB.ndiagL, AsB.ndiagU,
        AsB.Ainv, AsB.lda, AsB.ipiv, AsB.info);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (ARTYPE shift).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::
FactorAsB(ARFLOAT sigmaR, ARFLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Defining matrix AsBc.

  if (!AsBc.IsDefined()) {
    part = partp;
    AsBc.DefineMatrix(A->ncols(), max(A->ndiagL,B->ndiagL),
                      max(A->ndiagU,B->ndiagU), 0);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsBc.CreateStructure();

  // Subtracting sigma*B from A and storing the result on AsBc.

  SubtractAsB(sigmaR, sigmaI);

  // Decomposing AsBc.

  gbtrf(AsBc.n, AsBc.n, AsBc.ndiagL, AsBc.ndiagU, 
        AsBc.Ainv, AsBc.lda, AsBc.ipiv, AsBc.info);

  // Handling errors.

  AsBc.ThrowError();

  AsBc.factored = true;

} // FactorAsB (arcomplex<ARFLOAT> shift).
#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::
MultInvAsBv(arcomplex<ARFLOAT>* v, arcomplex<ARFLOAT>* w)
{

  AsB.MultInvv((ARTYPE*)v, (ARTYPE*)w);

} // MultInvAsBv (arcomplex<ARFLOAT>).
#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymPencil<ARTYPE, ARFLOAT>::MultInvAsBv(ARFLOAT* v, ARFLOAT* w)
{

  if (part == 'N') {    // shift is real.

    AsB.MultInvv((ARTYPE*)v, (ARTYPE*)w);

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
inline void ARbdNonSymPencil<ARTYPE, ARFLOAT>::
DefineMatrices(ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
               ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  A = &Ap;
  B = &Bp;

} // DefineMatrices.


template<class ARTYPE, class ARFLOAT>
inline ARbdNonSymPencil<ARTYPE, ARFLOAT>::
ARbdNonSymPencil(ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                 ARbdNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE, class ARFLOAT>
ARbdNonSymPencil<ARTYPE, ARFLOAT>& ARbdNonSymPencil<ARTYPE, ARFLOAT>::
operator=(const ARbdNonSymPencil<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBNSPEN_H
