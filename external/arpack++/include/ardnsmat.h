/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDNSMat.h.
   Arpack++ class ARdsNonSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "ardnspen.h"

#ifndef ARDNSMAT_H
#define ARDNSMAT_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "ardfmat.h"

template<class AR_T, class AR_S> class ARdsNonSymPencil;

template<class ARTYPE, class ARFLOAT>
class ARdsNonSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARdsNonSymPencil<ARTYPE, ARFLOAT>;
  friend class ARdsNonSymPencil<ARFLOAT, ARFLOAT>;

 protected:

  bool                factored;
  int                 info;
  int*                ipiv;
  ARTYPE*             A;
  ARTYPE*             Ainv;
  ARdfMatrix<ARTYPE>  mat;

  void ClearMem(); 

  virtual void Copy(const ARdsNonSymMatrix& other);

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

  void DefineMatrix(int np, ARTYPE* Ap);

  void DefineMatrix(int mp, int np, ARTYPE* Ap);

  ARdsNonSymMatrix(): ARMatrix<ARTYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARdsNonSymMatrix(int np, ARTYPE* Ap);
  // Long constructor (square matrix).

  ARdsNonSymMatrix(int mp, int np, ARTYPE* Ap);
  // Long constructor (rectangular matrix).

  ARdsNonSymMatrix(const std::string& file, int blksizep = 0);
  // Long constructor (Matrix stored in a file).

  ARdsNonSymMatrix(const ARdsNonSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsNonSymMatrix() { ClearMem(); }
  // Destructor.

  ARdsNonSymMatrix& operator=(const ARdsNonSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARdsNonSymMatrix member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
inline void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class ARTYPE, class ARFLOAT>
inline void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::
Copy(const ARdsNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  this->m         = other.m;
  this->n         = other.n;
  this->defined   = other.defined;
  factored  = other.factored;
  info      = other.info;
  A         = other.A;

  // Copying mat.

  if (other.mat.IsDefined()) {
    mat.Define(other.mat.Filename(),other.mat.BlockSize());
  }

  // Returning from here if "other" was not factored.

  if (!factored) return;

  // Copying vectors.

  Ainv = new ARTYPE[this->m*this->n];
  ipiv = new int[this->n];

  copy(this->m*this->n, other.Ainv, 1, Ainv, 1);
  for (int i=0; i<this->n; i++) ipiv[i] = other.ipiv[i];

} // Copy.


template<class ARTYPE, class ARFLOAT>
inline void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::CreateStructure()
{

  ClearMem();
  Ainv = new ARTYPE[this->m*this->n];
  ipiv = new int[this->n];

} // CreateStructure.


template<class ARTYPE, class ARFLOAT>
inline void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARdsNonSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARdsNonSymMatrix::FactorA");
  }

} // ThrowError.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::FactorA()
{

  // Quitting the function if A was not defined or is rectangular.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsNonSymMatrix::FactorA");
  }

  if (this->m!=this->n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymMatrix::FactorA");
  }

  if (mat.IsOutOfCore()) {
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARdsNonSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ::copy(this->m*this->n, A, 1, Ainv, 1);

  // Decomposing A.

  getrf(this->m, this->n, Ainv, this->m, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined or is rectangular.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARdsNonSymMatrix::FactorAsI");
  }

  if (this->m!=this->n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymMatrix::FactorAsI");
  }

  if (mat.IsOutOfCore()) {
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARdsNonSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  ::copy(this->m*this->n,A,1,Ainv,1);
  for (int i=0; i<(this->m*this->n); i+=this->m+1) Ainv[i]-=sigma;

  // Decomposing AsI.

  getrf(this->m, this->n, Ainv, this->m, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::MultMv(ARTYPE* v, ARTYPE* w)
{

  int     i;
  ARTYPE* t;
  ARTYPE  one;
  ARTYPE  zero;

  one  = (ARTYPE)0 + 1.0;
  zero = (ARTYPE)0;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsNonSymMatrix::MultMv");
  }

  // Determining w = M.v.

  if (mat.IsOutOfCore()) {

    if (this->m>this->n) { 

      // Matrix is "tall".

      mat.Rewind();
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("N", mat.RowsInMemory(), this->n, one, mat.Entries(), 
             mat.RowsInMemory(), v, 1, zero, &w[mat.FirstIndex()], 1);
      }

    }
    else {

      // Matrix is "fat".

      mat.Rewind();
      t = new ARTYPE[mat.ColsInMemory()];
      for (i=0; i<this->m; i++) w[i] = zero;
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("N", this->m, mat.ColsInMemory(), one, mat.Entries(), 
             this->m, &v[mat.FirstIndex()], 1, zero, t, 1);
        axpy(this->m, one, t, 1, w, 1); 
      }
      delete[] t;

    }

  }
  else {

    gemv("N", this->m, this->n, one, A, this->m, v, 1, zero, w, 1);

  }

} // MultMv.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::MultMtv(ARTYPE* v, ARTYPE* w)
{

  int     i;
  ARTYPE* t;
  ARTYPE  one;   
  ARTYPE  zero; 

  one  = (ARTYPE)0 + 1.0;
  zero = (ARTYPE)0;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsNonSymMatrix::MultMtv");
  }

  // Determining w = M'.v.

  if (mat.IsOutOfCore()) {

    if (this->m<=this->n) { 

      // Matrix is "fat".

      mat.Rewind();
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("T", this->m, mat.ColsInMemory(), one, mat.Entries(), 
             this->m, v, 1, zero, &w[mat.FirstIndex()], 1);
      }

    }
    else {

      // Matrix is "tall".

      mat.Rewind();
      t = new ARTYPE[mat.ColsInMemory()];
      for (i=0; i<this->m; i++) w[i] = zero;
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("T", mat.RowsInMemory(), this->n, one, mat.Entries(), 
             mat.RowsInMemory(), &v[mat.FirstIndex()], 1, zero, t, 1);
        axpy(mat.RowsInMemory(), one, t, 1, w, 1); 
      }
      delete[] t;

    }

  }
  else {

    gemv("T", this->m, this->n, one, A, this->m, v, 1, zero, w, 1);

  }


} // MultMtv.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::MultMtMv(ARTYPE* v, ARTYPE* w)
{

  int    i;
  ARTYPE *t, *s;
  ARTYPE one;   
  ARTYPE zero; 

  one  = (ARTYPE)0 + 1.0;
  zero = (ARTYPE)0;

  if (mat.IsOutOfCore() && (this->m>this->n)) {

    // Special code for "tall" matrices.

    t = new ARTYPE[mat.BlockSize()];
    s = new ARTYPE[this->n];

    mat.Rewind();
    for (i=0; i<this->n; i++) w[i] = zero;
    for (i=0; i<mat.NBlocks(); i++) {
      mat.ReadBlock();
      gemv("N", mat.RowsInMemory(), this->n, one, mat.Entries(), 
           mat.RowsInMemory(), v, 1, zero, t, 1);
      gemv("T", mat.RowsInMemory(), this->n, one, mat.Entries(), 
           mat.RowsInMemory(), t, 1, zero, s, 1);
      axpy(this->n, one, s, 1, w, 1); 

    }

    delete[] t;
    delete[] s;

  }
  else {

    t = new ARTYPE[this->m];

    MultMv(v,t);
    MultMtv(t,w);

    delete[] t;

  }


} // MultMtMv.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::MultMMtv(ARTYPE* v, ARTYPE* w)
{

  int    i;
  ARTYPE *t, *s;
  ARTYPE one;   
  ARTYPE zero; 

  one  = (ARTYPE)0 + 1.0;
  zero = (ARTYPE)0;

  if (mat.IsOutOfCore() && (this->m<=this->n)) {

    // Special code for "fat" matrices.

    t = new ARTYPE[mat.BlockSize()];
    s = new ARTYPE[this->m];

    mat.Rewind();
    for (i=0; i<this->m; i++) w[i] = zero;
    for (i=0; i<mat.NBlocks(); i++) {
      mat.ReadBlock();
      gemv("T", this->m, mat.ColsInMemory(), one, mat.Entries(), 
           this->m, v, 1, zero, t, 1);
      gemv("N", this->m, mat.ColsInMemory(), one, mat.Entries(), 
           this->m, t, 1, zero, s, 1);
      axpy(this->m, one, s, 1, w, 1); 

    }

    delete[] t;
    delete[] s;

  }
  else {

    t = new ARTYPE[this->n];

    MultMtv(v,t);
    MultMv(t,w);

    delete[] t;

  }

} // MultMMtv.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::Mult0MMt0v(ARTYPE* v, ARTYPE* w)
{

  MultMv(&v[this->m],w);
  MultMtv(v,&w[this->m]);

} // Mult0MMt0v.


template<class ARTYPE, class ARFLOAT>
void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARdsNonSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy(this->n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  getrs("N", this->n, 1, Ainv, this->m, ipiv, w, this->m, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class ARTYPE, class ARFLOAT>
inline void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int np, ARTYPE* Ap)
{

  // Defining member variables.

  this->n         = np;
  this->m         = np;
  A         = Ap;
  this->defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix (square).


template<class ARTYPE, class ARFLOAT>
inline void ARdsNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int mp, int np, ARTYPE* Ap)
{

  // Defining member variables.

  this->m         = mp;
  this->n         = np;
  A         = Ap;
  this->defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0;

} // DefineMatrix (rectangular).


template<class ARTYPE, class ARFLOAT>
inline ARdsNonSymMatrix<ARTYPE, ARFLOAT>::
ARdsNonSymMatrix(int np, ARTYPE* Ap) : ARMatrix<ARTYPE>(np)
{

  factored = false;
  DefineMatrix(np, Ap);

} // Long constructor (square matrix).


template<class ARTYPE, class ARFLOAT>
inline ARdsNonSymMatrix<ARTYPE, ARFLOAT>::
ARdsNonSymMatrix(int mp, int np, ARTYPE* Ap) : ARMatrix<ARTYPE>(mp, np)
{

  factored = false;
  DefineMatrix(mp, np, Ap);

} // Long constructor (rectangular matrix).


template<class ARTYPE, class ARFLOAT>
ARdsNonSymMatrix<ARTYPE, ARFLOAT>::ARdsNonSymMatrix(const std::string& file, int blksizep)
{

  factored = false;

  try {
    mat.Define(file, blksizep);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARdsNonSymMatrix");
  }

  if (mat.NCols() == mat.NRows()) {
    DefineMatrix(mat.NCols(), (ARTYPE*)mat.Entries());
  }
  else {                             
    DefineMatrix(mat.NRows(), mat.NCols(), (ARTYPE*)mat.Entries());
  }

} // Long constructor (Matrix stored in a file).


template<class ARTYPE, class ARFLOAT>
ARdsNonSymMatrix<ARTYPE, ARFLOAT>& ARdsNonSymMatrix<ARTYPE, ARFLOAT>::
operator=(const ARdsNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDNSMAT_H
