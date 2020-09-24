/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDNSPen.h.
   Arpack++ class ARdsNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDNSPEN_H
#define ARDNSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "ardnsmat.h"


template<class ARTYPE, class ARFLOAT>
class ARdsNonSymPencil
{

 protected:

  char                               part;
  ARdsNonSymMatrix<ARTYPE, ARFLOAT>* A;
  ARdsNonSymMatrix<ARTYPE, ARFLOAT>* B;
  ARdsNonSymMatrix<ARTYPE, ARFLOAT>  AsB;
#ifdef ARCOMP_H
  ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> AsBc;
#endif

  virtual void Copy(const ARdsNonSymPencil& other);

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

  void DefineMatrices(ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                      ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Bp);

  ARdsNonSymPencil() { part = 'N'; }
  // Short constructor that does nothing.

  ARdsNonSymPencil(ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                   ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Bp);
  // Long constructor.

  ARdsNonSymPencil(const ARdsNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsNonSymPencil() { }
  // Destructor.

  ARdsNonSymPencil& operator=(const ARdsNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARdsNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
inline void ARdsNonSymPencil<ARTYPE, ARFLOAT>::
Copy(const ARdsNonSymPencil<ARTYPE, ARFLOAT>& other)
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
void ARdsNonSymPencil<ARTYPE, ARFLOAT>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Copying A to AsB if sigma = 0.

  if (sigma == (ARTYPE)0) {

    AsB = *A;
    if (!AsB.IsFactored()) AsB.FactorA();
    return;

  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {
    AsB.DefineMatrix(A->ncols(), A->A);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure();

  // Subtracting sigma*B from A and storing the result on AsB.

  ::copy(A->m*A->n, A->A, 1, AsB.Ainv, 1); 
  axpy(A->m*A->n, -sigma, B->A, 1, AsB.Ainv, 1);

  // Decomposing AsB.

  getrf(AsB.m, AsB.n, AsB.Ainv, AsB.m, AsB.ipiv, AsB.info);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (ARTYPE shift).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARdsNonSymPencil<ARTYPE, ARFLOAT>::
FactorAsB(ARFLOAT sigmaR, ARFLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsBc.IsDefined()) {
    part = partp;
    AsBc.DefineMatrix(A->ncols(), 0);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsBc.CreateStructure(); 

  // Subtracting sigma*B from A and storing the result on AsBc.

  arcomplex<ARFLOAT> sigma(sigmaR, sigmaI);
  for (int i=0; i<(A->m*A->n); i++) AsBc.Ainv[i] = A->A[i]-sigma*B->A[i];

  // Decomposing AsBc.

  getrf(AsBc.m, AsBc.n, AsBc.Ainv, AsBc.m, AsBc.ipiv, AsBc.info);

  // Handling errors.

  AsBc.ThrowError();

  AsBc.factored = true;

} // FactorAsB (arcomplex<ARFLOAT> shift).
#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymPencil<ARTYPE, ARFLOAT>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H

template<class ARTYPE, class ARFLOAT>
void ARdsNonSymPencil<ARTYPE, ARFLOAT>::
MultInvAsBv(arcomplex<ARFLOAT>* v, arcomplex<ARFLOAT>* w)
{

  AsB.MultInvv((ARTYPE*)v,(ARTYPE*)w);

} // MultInvAsBv (arcomplex<ARFLOAT>).

#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymPencil<ARTYPE, ARFLOAT>::MultInvAsBv(ARFLOAT* v, ARFLOAT* w)
{

  if (part == 'N') {    // shift is real.

    AsB.MultInvv((ARTYPE*)v,(ARTYPE*)w);

  }
  else {                // shift is complex.

#ifdef ARCOMP_H

    int              i;
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
inline void ARdsNonSymPencil<ARTYPE, ARFLOAT>::
DefineMatrices(ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
               ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  A = &Ap;
  B = &Bp;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARdsNonSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE, class ARFLOAT>
inline ARdsNonSymPencil<ARTYPE, ARFLOAT>::
ARdsNonSymPencil(ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                 ARdsNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE, class ARFLOAT>
ARdsNonSymPencil<ARTYPE, ARFLOAT>& ARdsNonSymPencil<ARTYPE, ARFLOAT>::
operator=(const ARdsNonSymPencil<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDNSPEN_H
