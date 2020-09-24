/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDSPen.h.
   Arpack++ class ARdsSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDSPEN_H
#define ARDSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "ardsmat.h"


template<class ARTYPE>
class ARdsSymPencil
{

 protected:

  ARdsSymMatrix<ARTYPE>* A;
  ARdsSymMatrix<ARTYPE>* B;
  ARdsSymMatrix<ARTYPE>  AsB;

  virtual void Copy(const ARdsSymPencil& other);

  void SubtractAsB(ARTYPE sigma);

 public:

  bool IsFactored() { return AsB.IsFactored(); }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  void MultInvAsBv(ARTYPE* v, ARTYPE* w) {  AsB.MultInvv(v,w); }

  void DefineMatrices(ARdsSymMatrix<ARTYPE>& Ap, ARdsSymMatrix<ARTYPE>& Bp);

  ARdsSymPencil() { AsB.factored = false; }
  // Short constructor that does nothing.

  ARdsSymPencil(ARdsSymMatrix<ARTYPE>& Ap, ARdsSymMatrix<ARTYPE>& Bp);
  // Long constructor.

  ARdsSymPencil(const ARdsSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsSymPencil() { }
  // Destructor.

  ARdsSymPencil& operator=(const ARdsSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARdsSymPencil member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARdsSymPencil<ARTYPE>::Copy(const ARdsSymPencil<ARTYPE>& other)
{

  A        = other.A;
  B        = other.B;
  AsB      = other.AsB;

} // Copy.


template<class ARTYPE>
void ARdsSymPencil<ARTYPE>::SubtractAsB(ARTYPE sigma)
{

  int sizeA, i, j, k, l;

  // Copying A into AsB.

  sizeA = (A->ncols()*A->ncols()+A->ncols())/2;
  ::copy(sizeA, A->A, 1, AsB.Ainv, 1);

  // Returning if sigma == 0.

  if (sigma == (ARTYPE)0) return;

  // Subtracting sigma*B.

  if (A->uplo == B->uplo) {

    axpy(sizeA, -sigma, B->A, 1, AsB.Ainv, 1); 

  }
  else if (A->uplo == 'L') { // B->uplo == 'U'

    j = 0;
    for (i=0; i<A->n; i++) {
      for (l=i+1, k=(l*l+l)/2-1; l<=A->n; k+=(l++)) {
        AsB.Ainv[j++]-=sigma*B->A[k];
      }
    }

  }
  else { // A->uplo == 'U'  &&  B->uplo == 'L'

    j = 0;
    for (i=0; i<A->n; i++) {
      for (l=i+1, k=(l*l+l)/2-1; l<=A->n; k+=(l++)) {
        AsB.Ainv[k]-=sigma*B->A[j++];
      }
    }

  }
    
} // SubtractAsB (ARTYPE shift).


template<class ARTYPE>
void ARdsSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARdsSymPencil::FactorAsB");
  }

  // Copying A to AsB if sigma = 0.

  if (sigma == (ARTYPE)0) {

    AsB = *A;
    if (!AsB.IsFactored()) AsB.FactorA();
    return;

  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {
    AsB.DefineMatrix(A->ncols(), A->A, A->uplo);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure(); 

  // Subtracting sigma*B from A and storing the result on AsB.

  SubtractAsB(sigma);

  // Decomposing AsB.

  sptrf(&AsB.uplo, AsB.n, AsB.Ainv, AsB.ipiv, AsB.info);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (ARTYPE shift).


template<class ARTYPE>
void ARdsSymPencil<ARTYPE>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  copy(A->ncols(), w, 1, v, 1);
  B->MultInvv(w, w);

} // MultInvBAv.


template<class ARTYPE>
inline void ARdsSymPencil<ARTYPE>::
DefineMatrices(ARdsSymMatrix<ARTYPE>& Ap, ARdsSymMatrix<ARTYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARdsNonSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE>
inline ARdsSymPencil<ARTYPE>::
ARdsSymPencil(ARdsSymMatrix<ARTYPE>& Ap, ARdsSymMatrix<ARTYPE>& Bp)
{

  AsB.factored  = false;
  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE>
ARdsSymPencil<ARTYPE>& ARdsSymPencil<ARTYPE>::
operator=(const ARdsSymPencil<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDSPEN_H
