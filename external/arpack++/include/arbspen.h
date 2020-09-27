/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBSPen.h.
   Arpack++ class ARbdSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBSPEN_H
#define ARBSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arbsmat.h"


template<class ARTYPE>
class ARbdSymPencil
{

 protected:

  ARbdSymMatrix<ARTYPE>* A;
  ARbdSymMatrix<ARTYPE>* B;
  ARbdSymMatrix<ARTYPE>  AsB;

  int max(int a, int b) { return (a>b)?a:b; }

  int min(int a, int b) { return (a<b)?a:b; }

  virtual void Copy(const ARbdSymPencil& other);

  void SubtractAsB(ARTYPE sigma);

 public:

  bool IsFactored() { return AsB.IsFactored(); }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  void MultInvAsBv(ARTYPE* v, ARTYPE* w) {  AsB.MultInvv(v,w); }

  void DefineMatrices(ARbdSymMatrix<ARTYPE>& Ap, ARbdSymMatrix<ARTYPE>& Bp);

  ARbdSymPencil() { AsB.factored = false; }
  // Short constructor that does nothing.

  ARbdSymPencil(ARbdSymMatrix<ARTYPE>& Ap, ARbdSymMatrix<ARTYPE>& Bp);
  // Long constructor.

  ARbdSymPencil(const ARbdSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdSymPencil() { }
  // Destructor.

  ARbdSymPencil& operator=(const ARbdSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARbdSymPencil member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARbdSymPencil<ARTYPE>::Copy(const ARbdSymPencil<ARTYPE>& other)
{

  A        = other.A;
  B        = other.B;
  AsB      = other.AsB;

} // Copy.


template<class ARTYPE>
void ARbdSymPencil<ARTYPE>::SubtractAsB(ARTYPE sigma)
{

  int    i, n, minD, ndA, ndB, ndAsB, lda;
  ARTYPE negsig;

  negsig = -sigma;
  n      = AsB.n;
  ndA    = A->nsdiag;
  ndB    = B->nsdiag;
  ndAsB  = AsB.nsdiag;
  lda    = AsB.lda;

  // Expanding A.

  if (A->uplo == 'U') {

    // Copying the main diagonal of A.

    copy(n, &A->A[ndA], ndA+1, &AsB.Ainv[2*ndAsB], lda);

    // Copying the superdiagonals of A.

    for (i = 0; i < ndA; i++) {
      copy(n, &A->A[i], ndA+1, &AsB.Ainv[2*ndAsB-ndA+i], lda);
      copy(n-ndA+i, &A->A[i+(ndA-i)*(ndA+1)], ndA+1, 
           &AsB.Ainv[2*ndAsB+ndA-i], lda);
    }

  }
  else {

    // Copying the main diagonal of A to Ainv.

    copy(n, &A->A[0], ndA+1, &AsB.Ainv[2*ndAsB], lda);

    // Copying the subdiagonals of A to Ainv.

    for (i = 1; i <= ndA; i++) {
      copy(n, &A->A[i], ndA+1, &AsB.Ainv[2*ndAsB+i], lda);
      copy(n-i, &A->A[i], ndA+1, &AsB.Ainv[2*ndAsB-i+i*lda], lda);
    }

  }

  // Transferring part of B (*(-sigma)) if AsB.nsdiag > A->nsdiag.

  if (A->nsdiag < AsB.nsdiag) {

    if (B->uplo == 'U') {

      for (i = 0; i < ndAsB-ndA; i++) {
        copy(n, &B->A[i], ndB+1, &AsB.Ainv[ndAsB+i], lda);
        scal(n, negsig, &AsB.Ainv[ndAsB+i], lda);
        copy(n-ndAsB+i, &AsB.Ainv[ndAsB+i+(ndAsB-i)*lda], lda,
             &AsB.Ainv[lda-i-1], lda);  
      }

    }
    else {

      for (i = ndA+1; i <= ndAsB; i++) {
        copy(n, &B->A[i], ndB+1, &AsB.Ainv[2*ndAsB+i], lda);
        scal(n, negsig, &AsB.Ainv[2*ndAsB+i], lda);
        copy(n-i, &AsB.Ainv[2*ndAsB+i], lda, 
             &AsB.Ainv[2*ndAsB-i+i*lda], lda);
      }

    }

  } 

  // Subtracting sigma*B from A.

  minD   = min(ndA, ndB);

  if (B->uplo == 'U') {

    // Subtracting the main diagonal of B.

      axpy(n, negsig, &B->A[ndB], ndB+1, &AsB.Ainv[2*ndAsB], lda);

    // Subtracting the superdiagonals.

    for (i = 0; i < minD; i++) {
      axpy(n, negsig, &B->A[ndB-minD+i], ndB+1, 
           &AsB.Ainv[2*ndAsB-minD+i], lda);
      copy(n-minD+i, &AsB.Ainv[2*ndAsB-minD+i+(minD-i)*lda], lda,
           &AsB.Ainv[2*ndAsB+minD-i], lda);  
    }

  }
  else {

    // Subtracting the main diagonal of B.

      axpy(n, negsig, &B->A[0], ndB+1, &AsB.Ainv[2*ndAsB], lda);

    // Subtracting the subdiagonals.

    for (i = 1; i <= minD; i++) {
      axpy(n, negsig, &B->A[i], ndB+1, &AsB.Ainv[2*ndAsB+i], lda);
      copy(n-i, &AsB.Ainv[2*ndAsB+i], lda, 
           &AsB.Ainv[2*ndAsB-i+i*lda], lda);
    }
   
  }

} // SubtractAsB (ARTYPE shift).


template<class ARTYPE>
void ARbdSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARbdSymPencil::FactorAsB");
  }

  // Copying A to AsB if sigma = 0.

  if (sigma == (ARTYPE)0) {

    AsB = *A;
    if (!AsB.IsFactored()) AsB.FactorA();
    return;

  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {
    AsB.DefineMatrix(A->ncols(), max(A->nsdiag, B->nsdiag), A->A);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure(); 

  // Subtracting sigma*B from A and storing the result on AsB.

  SubtractAsB(sigma);

  // Decomposing AsB.

  gbtrf(AsB.n, AsB.n, AsB.nsdiag, AsB.nsdiag, 
        AsB.Ainv, AsB.lda, AsB.ipiv, AsB.info);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (ARTYPE shift).


template<class ARTYPE>
void ARbdSymPencil<ARTYPE>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  copy(A->ncols(), w, 1, v, 1);
  B->MultInvv(w, w);

} // MultInvBAv.


template<class ARTYPE>
inline void ARbdSymPencil<ARTYPE>::
DefineMatrices(ARbdSymMatrix<ARTYPE>& Ap, ARbdSymMatrix<ARTYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

} // DefineMatrices.


template<class ARTYPE>
inline ARbdSymPencil<ARTYPE>::
ARbdSymPencil(ARbdSymMatrix<ARTYPE>& Ap, ARbdSymMatrix<ARTYPE>& Bp)
{

  AsB.factored  = false;
  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE>
ARbdSymPencil<ARTYPE>& ARbdSymPencil<ARTYPE>::
operator=(const ARbdSymPencil<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBSPEN_H
