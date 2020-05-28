/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBSMat.h.
   Arpack++ class ARbdSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arbspen.h"

#ifndef ARBSMAT_H
#define ARBSMAT_H

#include <cstddef>
#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"

template<class ARTYPE> class ARbdSymPencil;

template<class ARTYPE>
class ARbdSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARbdSymPencil<ARTYPE>;

 protected:

  bool     factored;
  char     uplo;
  int      nsdiag;
  int      lda;
  int      info;
  int*     ipiv;
  ARTYPE*  A;
  ARTYPE*  Ainv;

  void ClearMem(); 

  virtual void Copy(const ARbdSymMatrix& other);

  void ExpandA();

  void SubtractAsI(ARTYPE sigma);

  void CreateStructure();

  void ThrowError();
  
 public:

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int nsdiagp, ARTYPE* Ap, char uplop = 'L');

  ARbdSymMatrix(): ARMatrix<ARTYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARbdSymMatrix(int np, int nsdiagp, ARTYPE* Ap, char uplop = 'L');
  // Long constructor.

  ARbdSymMatrix(const ARbdSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdSymMatrix() { ClearMem(); }
  // Destructor.

  ARbdSymMatrix& operator=(const ARbdSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARbdSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARbdSymMatrix<ARTYPE>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class ARTYPE>
inline void ARbdSymMatrix<ARTYPE>::
Copy(const ARbdSymMatrix<ARTYPE>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  this->m         = other.m;
  this->n         = other.n;
  this->defined   = other.defined;
  factored  = other.factored;
  uplo      = other.uplo;
  nsdiag    = other.nsdiag;
  lda       = other.lda;
  info      = other.info;
  A         = other.A;

  // Returning from here if "other" was not factored.

  if (!factored) return;

  // Copying vectors.

  Ainv = new ARTYPE[this->n*lda];
  ipiv = new int[this->n];

  copy(this->n*lda, other.Ainv, 1, Ainv, 1);
  for (int i=0; i<this->n; i++) ipiv[i] = other.ipiv[i];

} // Copy.


template<class ARTYPE>
void ARbdSymMatrix<ARTYPE>::ExpandA()
{

  int i;
 
  if (uplo == 'U') {

    // Copying the main diagonal of A to Ainv.

    copy(this->n, &A[nsdiag], nsdiag+1, &Ainv[2*nsdiag], lda);

    // Copying the superdiagonals of A to Ainv.

    for (i = 0; i < nsdiag; i++) {
      copy(this->n, &A[i], nsdiag+1, &Ainv[nsdiag+i], lda);
      copy(this->n-nsdiag+i, &A[i+(nsdiag-i)*(nsdiag+1)], nsdiag+1, 
           &Ainv[3*nsdiag-i], lda);
    }

  }
  else {

    // Copying the main diagonal of A to Ainv.

    copy(this->n, &A[0], nsdiag+1, &Ainv[2*nsdiag], lda);

    // Copying the subdiagonals of A to Ainv.

    for (i = 1; i <= nsdiag; i++) {
      copy(this->n, &A[i], nsdiag+1, &Ainv[2*nsdiag+i], lda);
      copy(this->n-i, &A[i], nsdiag+1, &Ainv[2*nsdiag-i+i*lda], lda);
    }

  }

} // ExpandA.


template<class ARTYPE>
void ARbdSymMatrix<ARTYPE>::SubtractAsI(ARTYPE sigma)
{

  // Copying A to Ainv.

  ExpandA();

  // Subtracting sigma from diagonal elements.

  for (int i=(2*nsdiag); i<(lda*this->n); i+=lda) Ainv[i] -= sigma; 

} // SubtractAsI.


template<class ARTYPE>
inline void ARbdSymMatrix<ARTYPE>::CreateStructure()
{

  ClearMem();
  Ainv = new ARTYPE[lda*this->n];
  ipiv = new int[this->n];

} // CreateStructure.


template<class ARTYPE>
inline void ARbdSymMatrix<ARTYPE>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARbdSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARbdSymMatrix::FactorA");
  }

} // ThrowError.


template<class ARTYPE>
void ARbdSymMatrix<ARTYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ExpandA();

  // Decomposing A.

  gbtrf(this->n, this->n, nsdiag, nsdiag, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class ARTYPE>
void ARbdSymMatrix<ARTYPE>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  gbtrf(this->n, this->n, nsdiag, nsdiag, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class ARTYPE>
void ARbdSymMatrix<ARTYPE>::MultMv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE  one  = (ARTYPE)0 + 1.0;
  ARTYPE  zero = (ARTYPE)0;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdSymMatrix::MultMv");
  }

  // Determining w = M.v.

  sbmv(&uplo, this->n, nsdiag, one, A, nsdiag+1, v, 1, zero, w, 1);

} // MultMv.


template<class ARTYPE>
void ARbdSymMatrix<ARTYPE>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARbdSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy(this->n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  gbtrs("N", this->n, nsdiag, nsdiag, 1, Ainv, lda, ipiv, w, this->m, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class ARTYPE>
inline void ARbdSymMatrix<ARTYPE>::
DefineMatrix(int np, int nsdiagp, ARTYPE* Ap, char uplop)
{

  // Defining member variables.

  this->m         = np;
  this->n         = np;
  nsdiag    = nsdiagp;
  lda       = 3*nsdiag+1;
  uplo      = uplop;
  A         = Ap;
  this->defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix.


template<class ARTYPE>
inline ARbdSymMatrix<ARTYPE>::
ARbdSymMatrix(int np, int nsdiagp, 
              ARTYPE* Ap, char uplop) : ARMatrix<ARTYPE>(np)
{

  factored = false;
  DefineMatrix(np, nsdiagp, Ap, uplop);

} // Long constructor.


template<class ARTYPE>
ARbdSymMatrix<ARTYPE>& ARbdSymMatrix<ARTYPE>::
operator=(const ARbdSymMatrix<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBSMAT_H
