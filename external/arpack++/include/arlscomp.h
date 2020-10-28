/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLSComp.h.
   Arpack++ class ARluCompStdEig definition
   (superlu version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLSCOMP_H
#define ARLSCOMP_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arscomp.h"
#include "arlnsmat.h"
#include "arrseig.h"


template<class ARFLOAT>
class ARluCompStdEig:
  public virtual ARCompStdEig<ARFLOAT, 
                              ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> > {

 protected:

 // a) Protected function:

  virtual void Copy(const ARluCompStdEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // b) Public functions:

 // b.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<ARFLOAT> sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<ARFLOAT> sigmap);

 // b.2) Constructors and destructor.

  ARluCompStdEig() { }
  // Short constructor.

  ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 const std::string& whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 arcomplex<ARFLOAT> sigma, const std::string& whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluCompStdEig(const ARluCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluCompStdEig() { }
  // Destructor.


 // c) Operators.

  ARluCompStdEig& operator=(const ARluCompStdEig& other);
  // Assignment operator.

}; // class ARluCompStdEig.


// ------------------------------------------------------------------------ //
// ARluCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluCompStdEig<ARFLOAT>::
Copy(const ARluCompStdEig<ARFLOAT>& other)
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>, 
           ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> >::
    Copy(other);
  if (this->mode > 2) this->objOP->FactorAsI(this->sigmaR);

} // Copy.


template<class ARFLOAT>
inline void ARluCompStdEig<ARFLOAT>::ChangeShift(arcomplex<ARFLOAT> sigmap)
{

  this->objOP->FactorAsI(sigmap);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmap);

} // ChangeShift.


template<class ARFLOAT>
inline void ARluCompStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>, 
           ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetRegularMode(this->objOP, 
                   &ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluCompStdEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>, 
           ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP,
                       &ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARluCompStdEig<ARFLOAT>::
ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               const std::string& whichp, int ncvp, ARFLOAT tolp,
               int maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A,
                   &ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluCompStdEig<ARFLOAT>::
ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               arcomplex<ARFLOAT> sigmap, const std::string& whichp, int ncvp,
               ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A,
                   &ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARluCompStdEig<ARFLOAT>& ARluCompStdEig<ARFLOAT>::
operator=(const ARluCompStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLSCOMP_H
