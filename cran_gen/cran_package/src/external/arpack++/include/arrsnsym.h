/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARRSNSym.h.
   Arpack++ class ARrcNonSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRSNSYM_H
#define ARRSNSYM_H

#include <cstddef>

#include "arch.h"
#include "arerror.h"
#include "debug.h"
#include "arrseig.h"
#include "naupp.h"
#include "neupp.h"


template<class ARFLOAT>
class ARrcNonSymStdEig: public virtual ARrcStdEig<ARFLOAT, ARFLOAT> {

 protected:

 // a) Protected functions:

 // a.1) Memory control functions.

  int ValSize() { return this->nev+1; }
  // Provides the size of array EigVal.

  void ValAllocate();
  // Creates arrays EigValR and EigValI.

  void WorkspaceAllocate();
  // Allocates workspace for nonsymmetric problems.


 // a.2) Functions that handle original FORTRAN ARPACK code.

  void Aupp();
  // Interface to FORTRAN subroutines SNAUPD and DNAUPD.

  void Eupp();
  // Interface to FORTRAN subroutines SNEUPD and DNEUPD.


 // a.3) Functions that check user defined parameters.

  int CheckNev(int nevp);
  // Does Range checking on nev.


 // a.4) Auxiliary functions required when using STL vector class.

  bool ConjEigVec(int i);
  // Indicates if EigVec[i] is the second eigenvector in 
  // a complex conjugate pair.

#ifdef ARCOMP_H
#ifdef STL_VECTOR_H

  vector<arcomplex<ARFLOAT> >* GenComplex(vector<ARFLOAT>* RealPart, 
                                          vector<ARFLOAT>* ImagPart, 
                                          bool conj = false);
  // Generates a complex vector Complex = RealPart + I*ImagPart
  // (or Complex = RealPart - I*ImagPart, if conj = true).

  vector<arcomplex<ARFLOAT> >* GenComplex(int dim, ARFLOAT* RealPart, 
                                          ARFLOAT* ImagPart, 
                                          bool conj = false);
  // Generates a complex vector Complex = RealPart + I*ImagPart
  // (or Complex = RealPart - I*ImagPart, if conj = true). dim
  // is the length of RealPart and ImagPart.

  vector<arcomplex<ARFLOAT> >* GenComplex(int dim, ARFLOAT* RealPart);
  // Generates a complex vector from a real vector. dim is the
  // length of RealPart.

#endif // STL_VECTOR_H.
#endif // ARCOMP_H.

 public:

 // b) Public functions:

 // b.1) Trace functions.

  void Trace(const int digit = -5, const int getv0 = 0, const int aupd = 1,
             const int aup2 = 0,  const int aitr = 0,  const int eigt = 0,
             const int apps = 0,  const int gets = 0,  const int eupd = 0)
  {
    nTraceOn(digit, getv0, aupd, aup2, aitr, eigt, apps, gets, eupd); 
  }
  // Turns on trace mode. 


 // b.2) Functions that permit step by step execution of ARPACK.

  ARFLOAT* GetVectorImag();
  // When ido = 3, this function indicates where the imaginary part
  // of the eigenvalues of the current Hessenberg matrix are located.


 // b.3) Functions that perform all calculations in one step.

  int Eigenvalues(ARFLOAT* &EigValRp, ARFLOAT* &EigValIp,
                  bool ivec = false, bool ischur = false);
  // Overrides arrays EigValRp with the real part and EigValIp 
  // with the imaginary part of the eigenvalues of the problem. 
  // Calculates eigenvectors and Schur vectors if requested.

  int EigenValVectors(ARFLOAT* &EigVecp, ARFLOAT* &EigValRp, 
                      ARFLOAT* &EigValIp, bool ischur = false);
  // Overrides array EigVecp sequentially with the eigenvectors of the
  // given eigen-problem. Also stores the eigenvalues in EigValRp and
  // EigValIp. Calculates Schur vectors if requested.


 // b.4) Functions that return elements of vectors and matrices.

#ifdef ARCOMP_H
  arcomplex<ARFLOAT> Eigenvalue(int i);
  // Furnishes i-eth eigenvalue.
#endif // ARCOMP_H.

  ARFLOAT EigenvalueReal(int i);
  // Provides the real part of the i-eth eigenvalue.

  ARFLOAT EigenvalueImag(int i);
  // Provides the imaginary part of the i-eth eigenvalue.

#ifdef ARCOMP_H
  arcomplex<ARFLOAT> Eigenvector(int i, int j);
  // Furnishes element j of the i-eth eigenvector.
#endif // ARCOMP_H.

  ARFLOAT EigenvectorReal(int i, int j);
  // Provides the real part of element j of the i-eth eigenvector.

  ARFLOAT EigenvectorImag(int i, int j);
  // Provides the imaginary part of element j of the i-eth eigenvector.


 // b.5) Functions that provide raw access to internal vectors and matrices.

  ARFLOAT* RawEigenvaluesImag();
  // Provides raw access to the imaginary part of eigenvalues.


 // b.6) Functions that use STL vector class.

#ifdef STL_VECTOR_H

#ifdef ARCOMP_H
  vector<arcomplex<ARFLOAT> >* StlEigenvalues(bool ivec = false, 
                                              bool ischur = false);
  // Calculates the eigenvalues and stores them in a single STL vector.
  // Also calculates eigenvectors and Schur vectors if requested.
#endif // ARCOMP_H.

  vector<ARFLOAT>* StlEigenvaluesReal();
  // Returns the real part of the eigenvalues.

  vector<ARFLOAT>* StlEigenvaluesImag();
  // Returns the imaginary part of the eigenvalues.

#ifdef ARCOMP_H
  vector<arcomplex<ARFLOAT> >* StlEigenvector(int i);
  // Returns the i-th eigenvector.
#endif // ARCOMP_H.

  vector<ARFLOAT>* StlEigenvectorReal(int i);
  // Returns the real part of the i-th eigenvector.

  vector<ARFLOAT>* StlEigenvectorImag(int i);
  // Returns the imaginary part of the i-th eigenvector.

#endif // STL_VECTOR_H.


 // b.7) Constructors and destructor.

  ARrcNonSymStdEig() { }
  // Short constructor.

  ARrcNonSymStdEig(int np, int nevp, const std::string& whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0, ARFLOAT* residp = NULL,
                   bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcNonSymStdEig(int np, int nevp, ARFLOAT sigma, const std::string& whichp = "LM",
                   int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcNonSymStdEig(const ARrcNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcNonSymStdEig() { }
  // Destructor.

 // c) Operators.

  ARrcNonSymStdEig& operator=(const ARrcNonSymStdEig& other);
  // Assignment operator.

}; // class ARrcNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARrcNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARrcNonSymStdEig<ARFLOAT>::ValAllocate()
{

  if (this->EigValR == NULL) {
    this->EigValR = new ARFLOAT[ValSize()];
    this->EigValI = new ARFLOAT[ValSize()];
    this->newVal = true;
  }

} // ValAllocate.


template<class ARFLOAT>
inline void ARrcNonSymStdEig<ARFLOAT>::WorkspaceAllocate()
{

  this->lworkl  = 3*this->ncv*(this->ncv+2);
  this->lworkv  = 3*this->ncv;
  this->lrwork  = 0;
  this->workl   = new ARFLOAT[this->lworkl+1];
  this->workv   = new ARFLOAT[this->lworkv+1];

} // WorkspaceAllocate.


template<class ARFLOAT>
inline void ARrcNonSymStdEig<ARFLOAT>::Aupp()
{

  naupp(this->ido,this-> bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n,
        this->iparam, this->ipntr, this->workd, this->workl, this->lworkl, this->info);

} // Aupp.


template<class ARFLOAT>
inline void ARrcNonSymStdEig<ARFLOAT>::Eupp()
{

  neupp(this->rvec, this->HowMny, this->EigValR, this->EigValI, this->EigVec, this->n, this->sigmaR,
        this->sigmaI, this->workv, this->bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V,
        this->n, this->iparam, this->ipntr, this->workd, this->workl, this->lworkl, this->info);

} // Eupp.


template<class ARFLOAT>
inline int ARrcNonSymStdEig<ARFLOAT>::CheckNev(int nevp)
{
    std::cout << nevp << "ddddd22g\n";fflush(stdout);
  if ((nevp<=1)||(nevp>=(this->n-1))) { // nev must satisfy 1 < nev < n-1.
      std::cout << nevp << "dddd21d22g\n";fflush(stdout);
    throw ArpackError(ArpackError::NEV_OUT_OF_BOUNDS);
  }
  return nevp;

} // CheckNev.


template<class ARFLOAT>
bool ARrcNonSymStdEig<ARFLOAT>::ConjEigVec(int i)
{

  if (this->EigValI[i] == (ARFLOAT)0.0) return false;
  int j = i-1;
  while ((j >= 0) && (this->EigValI[j] != (ARFLOAT)0.0)) j--;
  if (((i-j)%2) == 0) {
    return true;
  }
  else {
    return false;
  }

} // ConjEigVec.


#ifdef STL_VECTOR_H // Defining functions that use STL vector class.
#ifdef ARCOMP_H

template<class ARFLOAT>
vector<arcomplex<ARFLOAT> >* ARrcNonSymStdEig<ARFLOAT>::
GenComplex(vector<ARFLOAT>* RealPart, vector<ARFLOAT>* ImagPart, bool conj)
{

  // Defining variables.

  vector<arcomplex<ARFLOAT> >* Result;
  try {
    Result = new vector<arcomplex<ARFLOAT> >(ValSize());
  }
  catch (ArpackError) { return NULL; }
  ARFLOAT* rp  = RealPart->begin();
  ARFLOAT* ip  = ImagPart->begin();
  ARFLOAT* end = RealPart->end();
  arcomplex<ARFLOAT>* s = Result->begin();

  // Creating a complex vector.

  if (!conj) {
    while (rp != end) *s++ = arcomplex<ARFLOAT>(*rp++, *ip++);
  }
  else {
    while (rp != end) *s++ = arcomplex<ARFLOAT>(*rp++, -(*ip++));
  }

  return Result;

} // GenComplex (vector<ARFLOAT> version).


template<class ARFLOAT>
vector<arcomplex<ARFLOAT> >* ARrcNonSymStdEig<ARFLOAT>::
GenComplex(int dim, ARFLOAT* RealPart, ARFLOAT* ImagPart, bool conj)
{

  // Defining variables.

  vector<arcomplex<ARFLOAT> >* Result;
  try {
    Result = new vector<arcomplex<ARFLOAT> >(dim);
  }
  catch (ArpackError) { return NULL; }
  ARFLOAT* rp  = RealPart; 
  ARFLOAT* ip  = ImagPart; 
  ARFLOAT* end = &RealPart[dim];
  arcomplex<ARFLOAT>* s = Result->begin();

  // Creating a complex vector.

  if (!conj) {
    while (rp != end) *s++ = arcomplex<ARFLOAT>(*rp++, *ip++);
  }
  else {
    while (rp != end) *s++ = arcomplex<ARFLOAT>(*rp++, -(*ip++));
  }

  return Result;

} // GenComplex (ARFLOAT* version).


template<class ARFLOAT>
vector<arcomplex<ARFLOAT> >* ARrcNonSymStdEig<ARFLOAT>::
GenComplex(int dim, ARFLOAT* RealPart)
{

  // Defining variables.

  vector<arcomplex<ARFLOAT> >* Result;
  try {
    Result = new vector<arcomplex<ARFLOAT> >(dim);
  }
  catch (ArpackError) { return NULL; }
  ARFLOAT* rp  = RealPart; 
  ARFLOAT* end = &RealPart[dim];
  arcomplex<ARFLOAT>* s = Result->begin();

  // Copying a real vector into a complex vector.

  while (rp != end) *s++ = *rp++;

  return Result;

} // GenComplex.

#endif // ARCOMP_H.
#endif // STL_VECTOR_H.


template<class ARFLOAT>
ARFLOAT* ARrcNonSymStdEig<ARFLOAT>::GetVectorImag()
{

  if (this->ido != 3) {
    throw ArpackError(ArpackError::CANNOT_GET_VECTOR, "GetVectorImag");
  }
  return &this->workl[this->ipntr[6]];

} // GetVectorImag.


template<class ARFLOAT>
int ARrcNonSymStdEig<ARFLOAT>::
Eigenvalues(ARFLOAT* &EigValRp, ARFLOAT* &EigValIp, bool ivec, bool ischur)
{

  if (this->ValuesOK) {                                 // Eigenvalues are available.
    if ((EigValRp == NULL)&&(EigValIp == NULL)) { // Moving eigenvalues.
      EigValRp = this->EigValR;
      EigValIp = this->EigValI;
      this->EigValR  = NULL;
      this->EigValI  = NULL;
      this->newVal   = false;
      this->ValuesOK = false;
    }
    else {                                        // Copying eigenvalues.
      try {
        if (EigValRp == NULL) EigValRp = new ARFLOAT[ValSize()];
        if (EigValIp == NULL) EigValIp = new ARFLOAT[ValSize()];
      }
      catch (ArpackError) { return 0; }
      copy(this->nconv,this->EigValR,1,EigValRp,1);
      copy(this->nconv,this->EigValI,1,EigValIp,1);
    }
  }
  else {
    if (this->newVal) {
      delete[] this->EigValR;
      delete[] this->EigValI;
      this->newVal = false;
    }
    try {
      if (EigValRp == NULL) EigValRp = new ARFLOAT[ValSize()];
      if (EigValIp == NULL) EigValIp = new ARFLOAT[ValSize()];
    }
    catch (ArpackError) { return 0; }
    this->EigValR = EigValRp;
    this->EigValI = EigValIp;
    if (ivec) {                              // Finding eigenvalues and vectors.
      this->nconv = this->FindEigenvectors(ischur);
    }
    else {                                   // Finding eigenvalues only.
      this->nconv = this->FindEigenvalues();
    }
    this->EigValR = NULL;
    this->EigValI = NULL;
  }
  return this->nconv;

} // Eigenvalues(EigValRp, EigValIp, ivec, ischur).


template<class ARFLOAT>
int ARrcNonSymStdEig<ARFLOAT>::
EigenValVectors(ARFLOAT* &EigVecp, ARFLOAT* &EigValRp,
                ARFLOAT* &EigValIp, bool ischur)
{

  if (this->ValuesOK) {               // Eigenvalues are already available .
    this->nconv = Eigenvalues(EigValRp, EigValIp, false);
    this->nconv = this->Eigenvectors(EigVecp, ischur);
  }
  else {                        // Eigenvalues ans vectors are not available.
    if (this->newVec) {
      delete[] this->EigVec;
      this->newVec = false;
    }
    if (this->newVal) {
      delete[] this->EigValR;
      delete[] this->EigValI;
      this->newVal = false;
    }
    try {
      if (EigVecp  == NULL) EigVecp  = new ARFLOAT[ValSize()*this->n];
      if (EigValRp == NULL) EigValRp = new ARFLOAT[ValSize()];
      if (EigValIp == NULL) EigValIp = new ARFLOAT[ValSize()];
    }
    catch (ArpackError) { return 0; }
    this->EigVec  = EigVecp;
    this->EigValR = EigValRp;
    this->EigValI = EigValIp;
    this->nconv   = this->FindEigenvectors(ischur);
    this->EigVec  = NULL;
    this->EigValR = NULL;
    this->EigValI = NULL;
  }
  return this->nconv;

} // EigenValVectors(EigVecp, EigValRp, EigValIp, ischur).


#ifdef ARCOMP_H
template<class ARFLOAT>
inline arcomplex<ARFLOAT> ARrcNonSymStdEig<ARFLOAT>::Eigenvalue(int i)
{

  // Returning i-eth eigenvalue.

  if (!this->ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "Eigenvalue(i)");
  }
  else if ((i>=this->nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvalue(i)");
  }
  return arcomplex<ARFLOAT>(this->EigValR[i],this->EigValI[i]);

} // Eigenvalue(i).
#endif // ARCOMP_H


template<class ARFLOAT>
inline ARFLOAT ARrcNonSymStdEig<ARFLOAT>::EigenvalueReal(int i)
{

  // Returning the real part of i-eth eigenvalue.

  if (!this->ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "EigenvalueReal(i)");
  }
  else if ((i>=this->nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvalueReal(i)");
  }
  return this->EigValR[i];

} // EigenvalueReal(i).


template<class ARFLOAT>
inline ARFLOAT ARrcNonSymStdEig<ARFLOAT>::EigenvalueImag(int i)
{

  // Returning the imaginary part of i-eth eigenvalue.

  if (!this->ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "EigenvalueImag(i)");
  }
  else if ((i>=this->nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvalueImag(i)");
  }
  return this->EigValI[i];

} // EigenvalueImag(i).


#ifdef ARCOMP_H
template<class ARFLOAT>
inline arcomplex<ARFLOAT> ARrcNonSymStdEig<ARFLOAT>::
Eigenvector(int i, int j)
{

  // Returning element j of i-eth eigenvector.

  if ((!this->VectorsOK)||(!this->ValuesOK)) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "Eigenvector(i,j)");
  }
  else if ((i>=this->nconv)||(i<0)||(j>=this->n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvector(i,j)");
  }
  if (this->EigValI[i]==(ARFLOAT)0.0) {   // Real eigenvalue.
    return arcomplex<ARFLOAT>(this->EigVec[i*this->n+j],(ARFLOAT)0.0);
  }
  else {                            // Complex eigenvalue.
    if (this->EigValI[i]>(ARFLOAT)0.0) {  // with positive imaginary part.
      return arcomplex<ARFLOAT>(this->EigVec[i*this->n+j], this->EigVec[(i+1)*this->n+j]);
    }
    else {                          // with negative imaginary part.
      return arcomplex<ARFLOAT>(this->EigVec[(i-1)*this->n+j], -this->EigVec[i*this->n+j]);
    }
  }

} // Eigenvector(i,j).
#endif // ARCOMP_H


template<class ARFLOAT>
inline ARFLOAT ARrcNonSymStdEig<ARFLOAT>::EigenvectorReal(int i, int j)
{

  // Returning the real part of element j of i-eth eigenvector.

  if (!this->VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "EigenvectorReal(i,j)");
  }
  else if ((i>=this->nconv)||(i<0)||(j>=this->n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvectorReal(i,j)");
  }
  return this->EigVec[i*this->n+j];

} // EigenvectorReal(i,j).


template<class ARFLOAT>
inline ARFLOAT ARrcNonSymStdEig<ARFLOAT>::EigenvectorImag(int i, int j)
{

  // Returning the imaginary part of element j of i-eth eigenvector.

  if ((!this->VectorsOK)||(!this->ValuesOK)) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "EigenvectorImag(i,j)");
  }
  else if ((i>=this->nconv)||(i<0)||(j>=this->n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvectorImag(i,j)");
  }
  if (this->EigValI[i]==(ARFLOAT)0.0) {   // Real eigenvalue.
    return (ARFLOAT)0.0;
  }
  else {                            // Complex eigenvalue.
    if (this->EigValI[i]>(ARFLOAT)0.0) {  // with positive imaginary part.
      return this->EigVec[(i+1)*this->n+j];
    }
    else {                          // with negative imaginary part.
      return -this->EigVec[i*this->n+j];
    }
  }

} // EigenvectorImag(i,j).


template<class ARFLOAT>
inline ARFLOAT* ARrcNonSymStdEig<ARFLOAT>::RawEigenvaluesImag()
{

  if (!this->ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "RawEigenvaluesImag");
  }
  return this->EigValI;

} // RawEigenvaluesImag.


#ifdef STL_VECTOR_H // Defining some functions that use STL vector class.

#ifdef ARCOMP_H
template<class ARFLOAT>
inline vector<arcomplex<ARFLOAT> >* ARrcNonSymStdEig<ARFLOAT>::
StlEigenvalues(bool ivec, bool ischur)
{

  // Returning the eigenvalues in a STL vector.

  // Defining temporary variables.

  vector<ARFLOAT>* StlEigValR;
  vector<ARFLOAT>* StlEigValI;
  ARFLOAT*         ValRPtr;
  ARFLOAT*         ValIPtr;

  try {
    StlEigValR = new vector<ARFLOAT>(ValSize());
    StlEigValI = new vector<ARFLOAT>(ValSize());
  }
  catch (ArpackError) { return NULL; }

  // Finding Eigenvalues.

  ValRPtr = StlEigValR->begin();
  ValIPtr = StlEigValI->begin();
  nconv = Eigenvalues(ValRPtr, ValIPtr, ivec, ischur);
  vector<arcomplex<ARFLOAT> >* Val = GenComplex(StlEigValR, StlEigValI);

  // Deleting temporary variables.

  delete StlEigValR;
  delete StlEigValI;

  return Val;

} // StlEigenvalues.
#endif // ARCOMP_H.


template<class ARFLOAT>
inline vector<ARFLOAT>* ARrcNonSymStdEig<ARFLOAT>::StlEigenvaluesReal()
{

  // Returning the real part of the eigenvalues in a STL vector.

  vector<ARFLOAT>* StlEigValR;
  
  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "StlEigenvaluesReal");
  }
  try {
    StlEigValR = new vector<ARFLOAT>(EigValR, &EigValR[ValSize()]);
  }
  catch (ArpackError) { return NULL; }
  return StlEigValR;

} // StlEigenvaluesReal.


template<class ARFLOAT>
inline vector<ARFLOAT>* ARrcNonSymStdEig<ARFLOAT>::StlEigenvaluesImag()
{

  // Returning the imaginary part of the eigenvalues in a STL vector.

  vector<ARFLOAT>* StlEigValI;

  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "StlEigenvaluesImag");
  }
  try {
    StlEigValI = new vector<ARFLOAT>(EigValI, &EigValI[ValSize()]);
  }
  catch (ArpackError) { return NULL; }
  return StlEigValI;

} // StlEigenvaluesImag.


#ifdef ARCOMP_H
template<class ARFLOAT>
inline vector<arcomplex<ARFLOAT> >* ARrcNonSymStdEig<ARFLOAT>::
StlEigenvector(int i)
{

  // Returning the i-th eigenvector in a STL vector.

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvector(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvector(i)");
  }
  if (EigValI[i] == (ARFLOAT)0.0) { // Real eigenvector.
    return GenComplex(n, &EigVec[i*n]);
  }
  else if (!ConjEigVec(i)) {      // First eigenvector in a conjugate pair.
    return GenComplex(n, &EigVec[i*n], &EigVec[(i+1)*n]);
  }
  else {                          // Second eigenvector in a conjugate pair.
    return GenComplex(n, &EigVec[(i-1)*n], &EigVec[i*n], true);
  }

} // StlEigenvector(i).
#endif // ARCOMP_H.


template<class ARFLOAT>
inline vector<ARFLOAT>* ARrcNonSymStdEig<ARFLOAT>::StlEigenvectorReal(int i)
{

  // Returning the real part of the i-th eigenvector in a STL vector.

  vector<ARFLOAT>* Vec;

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvectorReal(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvectorReal(i)");
  }
  if (!ConjEigVec(i)) { // Real eigenvector or first in a conj. pair.
    try {
      Vec = new vector<ARFLOAT>(&EigVec[i*n], &EigVec[(i+1)*n]);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }
  else {                // Second eigenvector in a conjugate pair.
    try {
      Vec = new vector<ARFLOAT>(&EigVec[(i-1)*n], &EigVec[i*n]);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }

} // StlEigenvectorReal(i).


template<class ARFLOAT>
inline vector<ARFLOAT>* ARrcNonSymStdEig<ARFLOAT>::StlEigenvectorImag(int i)
{

  // Returning the imaginary part of the i-th eigenvector in a STL vector.

  vector<ARFLOAT>* Vec;

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvectorImag(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvectorImag(i)");
  }
  if (EigValI[i] == (ARFLOAT)0.0) { // Real eigenvector.
    try {
      Vec = new vector<ARFLOAT>(ValSize(), (ARFLOAT)0.0);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }
  else if (!ConjEigVec(i)) {      // First eigenvector in a conjugate pair.
    try {
      Vec = new vector<ARFLOAT>(&EigVec[(i+1)*n], &EigVec[(i+2)*n]);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }
  else {                          // Second eigenvector in a conjugate pair.
    try {
      Vec = new vector<ARFLOAT>(&EigVec[i*n], &EigVec[(i+1)*n]);
    }
    catch (ArpackError) { return NULL; }
    for (ARFLOAT* s = Vec->begin(); s != Vec->end(); s++) *s = -(*s);
    return Vec;
  }

} // StlEigenvectorImag(i).

#endif // STL_VECTOR_H.


template<class ARFLOAT>
inline ARrcNonSymStdEig<ARFLOAT>::
ARrcNonSymStdEig(int np, int nevp, const std::string& whichp, int ncvp,
                 ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARrcNonSymStdEig<ARFLOAT>::
ARrcNonSymStdEig(int np, int nevp, ARFLOAT sigmap, const std::string& whichp, int ncvp,
                 ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->ChangeShift(sigmap);
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARrcNonSymStdEig<ARFLOAT>& ARrcNonSymStdEig<ARFLOAT>::
operator=(const ARrcNonSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRSNSYM_H

