/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDSMat.h.
   Arpack++ class ARdsSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "ardspen.h"

#ifndef ARDSMAT_H
#define ARDSMAT_H

#include <cstddef>

#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"

template<class ARTYPE> class ARdsSymPencil;

template<class ARTYPE>
class ARdsSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARdsSymPencil<ARTYPE>;

 protected:

  bool     factored;
  char     uplo;
  int      info;
  int*     ipiv;
  ARTYPE*  A;
  ARTYPE*  Ainv;

  void ClearMem(); 

  virtual void Copy(const ARdsSymMatrix& other);

  void SubtractAsI(ARTYPE sigma);

  void CreateStructure();

  void ThrowError();
  
 public:

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, ARTYPE* Ap, char uplop = 'L');

  ARdsSymMatrix(): ARMatrix<ARTYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARdsSymMatrix(int np, ARTYPE* Ap, char uplop = 'L');
  // Long constructor.

  ARdsSymMatrix(const ARdsSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsSymMatrix() { ClearMem(); }
  // Destructor.

  ARdsSymMatrix& operator=(const ARdsSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARdsSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARdsSymMatrix<ARTYPE>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class ARTYPE>
inline void ARdsSymMatrix<ARTYPE>::
Copy(const ARdsSymMatrix<ARTYPE>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  this->m         = other.m;
  this->n         = other.n;
  this->defined   = other.defined;
  factored  = other.factored;
  uplo      = other.uplo;
  info      = other.info;
  A         = other.A;

  // Returning from here if "other" was not factored.

  if (!factored) return;

  // Copying vectors.

  Ainv = new ARTYPE[(this->n*this->n+this->n)/2];
  ipiv = new int[this->n];

  copy((this->n*this->n+this->n)/2, other.Ainv, 1, Ainv, 1);
  for (int i=0; i<this->n; i++) ipiv[i] = other.ipiv[i];

} // Copy.


template<class ARTYPE>
void ARdsSymMatrix<ARTYPE>::SubtractAsI(ARTYPE sigma)
{

  int i,j;

  // Copying A to Ainv.

  ::copy((this->n*this->n+this->n)/2 ,A, 1, Ainv, 1);

  // Subtracting sigma from diagonal elements.

  if (uplo=='L') {
    for (i=0, j=0; i<this->n; j+=(this->n-(i++))) Ainv[j] -= sigma;
  }
  else {
    for (i=0, j=0; i<this->n; j+=(++i)) Ainv[j] -= sigma;
  }

} // SubtractAsI.


template<class ARTYPE>
inline void ARdsSymMatrix<ARTYPE>::CreateStructure()
{

  ClearMem();
  Ainv = new ARTYPE[(this->n*this->n+this->n)/2];
  ipiv = new int[this->n];

} // CreateStructure.


template<class ARTYPE>
inline void ARdsSymMatrix<ARTYPE>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARdsSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARdsSymMatrix::FactorA");
  }

} // ThrowError.


template<class ARTYPE>
void ARdsSymMatrix<ARTYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ::copy((this->n*this->n+this->n)/2 ,A, 1, Ainv, 1);

  // Decomposing A.

  sptrf(&uplo, this->n, Ainv, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class ARTYPE>
void ARdsSymMatrix<ARTYPE>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  sptrf(&uplo, this->n, Ainv, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class ARTYPE>
void ARdsSymMatrix<ARTYPE>::MultMv(ARTYPE* v, ARTYPE* w)
{

  int     i, j;

  ARTYPE  zero = (ARTYPE)0;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsSymMatrix::MultMv");
  }

  // Determining w = M.v (unfortunately, the BLAS does not 
  // have a routine that works with packed matrices).

  for (i=0; i<this->n; i++) w[i] = zero;

  if (uplo=='L') {
 
    for (i=0, j=0; i<this->n; j+=(this->n-(i++))) {
      w[i] += dot(this->n-i, &A[j], 1, &v[i], 1);
      axpy(this->n-i-1, v[i], &A[j+1], 1, &w[i+1], 1);
    }
 
  }  
  else { // uplo = 'U'

    for (i=0, j=0; i<this->n; j+=(++i)) {
      w[i] += dot(i+1, &A[j], 1, v, 1);
      axpy(i, v[i], &A[j], 1, w, 1);
    }
   
  }

} // MultMv.


template<class ARTYPE>
void ARdsSymMatrix<ARTYPE>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARdsSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy(this->n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  sptrs(&uplo, this->n, 1, Ainv, ipiv, w, this->n, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class ARTYPE>
inline void ARdsSymMatrix<ARTYPE>::
DefineMatrix(int np, ARTYPE* Ap, char uplop)
{

  // Defining member variables.

  this->m         = np;
  this->n         = np;
  uplo      = uplop;
  A         = Ap;
  this->defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix.


template<class ARTYPE>
inline ARdsSymMatrix<ARTYPE>::
ARdsSymMatrix(int np, ARTYPE* Ap, char uplop) : ARMatrix<ARTYPE>(np)
{

  factored = false;
  DefineMatrix(np, Ap, uplop);

} // Long constructor.


template<class ARTYPE>
ARdsSymMatrix<ARTYPE>& ARdsSymMatrix<ARTYPE>::
operator=(const ARdsSymMatrix<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDSMAT_H
