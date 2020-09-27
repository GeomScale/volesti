/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBNSMat.h.
   Arpack++ class ARbdNonSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arbnspen.h"

#ifndef ARBNSMAT_H
#define ARBNSMAT_H

#include <cstddef>
#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"

template<class AR_T, class AR_S> class ARbdNonSymPencil;

template<class ARTYPE, class ARFLOAT>
class ARbdNonSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARbdNonSymPencil<ARTYPE, ARFLOAT>;
  friend class ARbdNonSymPencil<ARFLOAT, ARFLOAT>;

 protected:

  bool     factored;
  int      ndiagL;
  int      ndiagU;
  int      lda;
  int      info;
  int*     ipiv;
  ARTYPE*  A;
  ARTYPE*  Ainv;

  void ClearMem(); 

  virtual void Copy(const ARbdNonSymMatrix& other);

  void ExpandA();

  void SubtractAsI(ARTYPE sigma);

  void CreateStructure();

  void ThrowError();
  
 public:

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultMtv(ARTYPE* v, ARTYPE* w);

  void MultMtMv(ARTYPE* v, ARTYPE* w);

  void MultMMtv(ARTYPE* v, ARTYPE* w);

  void Mult0MMt0v(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int ndiagLp, int ndiagUp, ARTYPE* Ap);

  ARbdNonSymMatrix(): ARMatrix<ARTYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARbdNonSymMatrix(int np, int ndiagLp, int ndiagUp, ARTYPE* Ap);
  // Long constructor.

  ARbdNonSymMatrix(const ARbdNonSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdNonSymMatrix() { ClearMem(); }
  // Destructor.

  ARbdNonSymMatrix& operator=(const ARbdNonSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARbdNonSymMatrix member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
inline void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class ARTYPE, class ARFLOAT>
inline void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::
Copy(const ARbdNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  this->m         = other.m;
   this->n         = other.n;
   this->defined   = other.defined;
  factored  = other.factored;
  ndiagL    = other.ndiagL;
  ndiagU    = other.ndiagU;
  lda       = other.lda;
  info      = other.info;
  A         = other.A;

  // Returning from here if "other" was not factored.

  if (!factored) return;

  // Copying vectors.

  Ainv = new ARTYPE[ this->n*lda];
  ipiv = new int[ this->n];

  copy( this->n*lda, other.Ainv, 1, Ainv, 1);
  for (int i=0; i< this->n; i++) ipiv[i] = other.ipiv[i];

} // Copy.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::ExpandA()
{

  int i, inca;
 
  // Copying A to Ainv.

  inca = ndiagL+ndiagU+1;
  for (i = 0; i < inca; i++) {
    copy( this->n, &A[i], inca, &Ainv[ndiagL+i], lda);
  }

} // ExpandA.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::SubtractAsI(ARTYPE sigma)
{

  // Copying A to Ainv.

  ExpandA();

  // Subtracting sigma from diagonal elements.

  for (int i=(ndiagL+ndiagU); i<(lda* this->n); i+=lda) Ainv[i] -= sigma; 

} // SubtractAsI.


template<class ARTYPE, class ARFLOAT>
inline void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::CreateStructure()
{

  ClearMem();
  Ainv = new ARTYPE[lda* this->n];
  ipiv = new int[ this->n];

} // CreateStructure.


template<class ARTYPE, class ARFLOAT>
inline void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARbdNonSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARbdNonSymMatrix::FactorA");
  }

} // ThrowError.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::FactorA()
{

  // Quitting the function if A was not defined.

  if (! this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdNonSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ExpandA();

  // Decomposing A.

  gbtrf( this->n,  this->n, ndiagL, ndiagU, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.

  if (! this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARbdNonSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  gbtrf( this->n,  this->n, ndiagL, ndiagU, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::MultMv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE  one;
  ARTYPE  zero;

  one  = (ARTYPE)0 + 1.0;
  zero = (ARTYPE)0;

  // Quitting the function if A was not defined.

  if (! this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdNonSymMatrix::MultMv");
  }

  // Determining w = M.v.

  gbmv("N",  this->m,  this->n, ndiagL, ndiagU, one, A,
       ndiagL+ndiagU+1, v, 1, zero, w, 1);

} // MultMv.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::MultMtv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE  one;   
  ARTYPE  zero; 

  one  = (ARTYPE)0 + 1.0;
  zero = (ARTYPE)0;

  // Quitting the function if A was not defined.

  if (! this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdNonSymMatrix::MultMtv");
  }

  // Determining w = M'.v.

  gbmv("T",  this->m,  this->n, ndiagL, ndiagU, one, A,
       ndiagL+ndiagU+1, v, 1, zero, w, 1);   

} // MultMtv.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::MultMtMv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE* t = new ARTYPE[ this->m];

  MultMv(v,t);
  MultMtv(t,w);

  delete[] t;

} // MultMtMv.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::MultMMtv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE* t = new ARTYPE[ this->n];

  MultMtv(v,t);
  MultMv(t,w);

  delete[] t;

} // MultMMtv.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::Mult0MMt0v(ARTYPE* v, ARTYPE* w)
{

  MultMv(&v[ this->m],w);
  MultMtv(v,&w[ this->m]);

} // Mult0MMt0v.


template<class ARTYPE, class ARFLOAT>
void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARbdNonSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy( this->n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  gbtrs("N",  this->n, ndiagL, ndiagU, 1, Ainv, lda, ipiv, w,  this->m, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class ARTYPE, class ARFLOAT>
inline void ARbdNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int np, int ndiagLp, int ndiagUp, ARTYPE* Ap)
{

  // Defining member variables.

   this->m         = np;
   this->n         = np;
  ndiagL    = ndiagLp;
  ndiagU    = ndiagUp;
  lda       = 2*ndiagL+ndiagU+1;
  A         = Ap;
   this->defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix.


template<class ARTYPE, class ARFLOAT>
inline ARbdNonSymMatrix<ARTYPE, ARFLOAT>::
ARbdNonSymMatrix(int np, int ndiagLp, 
                 int ndiagUp, ARTYPE* Ap) : ARMatrix<ARTYPE>(np)
{

  factored = false;
  DefineMatrix(np, ndiagLp, ndiagUp, Ap);

} // Long constructor.


template<class ARTYPE, class ARFLOAT>
ARbdNonSymMatrix<ARTYPE, ARFLOAT>& ARbdNonSymMatrix<ARTYPE, ARFLOAT>::
operator=(const ARbdNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
     this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBNSMAT_H
