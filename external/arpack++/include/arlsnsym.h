/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLSNSym.h.
   Arpack++ class ARluNonSymStdEig definition
   (SuperLU version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLSNSYM_H
#define ARLSNSYM_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arsnsym.h"
#include "arlnsmat.h"

template<class ARFLOAT>
class ARluNonSymStdEig:
  public virtual ARNonSymStdEig<ARFLOAT, ARluNonSymMatrix<ARFLOAT, ARFLOAT> > {

 protected:

 // a) Protected function:

  virtual void Copy(const ARluNonSymStdEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // b) Public functions:

 // b.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // b.2) Constructors and destructor.

  ARluNonSymStdEig() { }
  // Short constructor.

  ARluNonSymStdEig(int nevp, ARluNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   const std::string& whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluNonSymStdEig(int nevp, ARluNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARFLOAT sigma, const std::string& whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluNonSymStdEig(const ARluNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluNonSymStdEig() { }
  // Destructor.

 // c) Operators.

  ARluNonSymStdEig& operator=(const ARluNonSymStdEig& other);
  // Assignment operator.

}; // class ARluNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARluNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluNonSymStdEig<ARFLOAT>::
Copy(const ARluNonSymStdEig<ARFLOAT>& other)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARluNonSymMatrix<ARFLOAT, ARFLOAT> >:: Copy(other);
  if (this->mode > 2) this->objOP->FactorAsI(this->sigmaR);

} // Copy.


template<class ARFLOAT>
inline void ARluNonSymStdEig<ARFLOAT>::ChangeShift(ARFLOAT sigmaRp)
{

  this->sigmaR    = sigmaRp;
  this->sigmaI    = 0.0;
  this->mode      = 3;
  this->iparam[7] = this->mode;

  this->objOP->FactorAsI(this->sigmaR);
  this->Restart();

} // ChangeShift.


template<class ARFLOAT>
inline void ARluNonSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARluNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetRegularMode(this->objOP, &ARluNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluNonSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARluNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP, 
                       &ARluNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARluNonSymStdEig<ARFLOAT>::
ARluNonSymStdEig(int nevp, ARluNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 const std::string& whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A, 
                   &ARluNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluNonSymStdEig<ARFLOAT>::
ARluNonSymStdEig(int nevp, ARluNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARFLOAT sigmap, const std::string& whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A, 
                   &ARluNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARluNonSymStdEig<ARFLOAT>& ARluNonSymStdEig<ARFLOAT>::
operator=(const ARluNonSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLSNSYM_H
