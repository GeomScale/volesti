/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLGComp.h.
   Arpack++ class ARluCompGenEig definition
   (superlu version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLGCOMP_H
#define ARLGCOMP_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arlnsmat.h"
#include "arlnspen.h"
#include "arrseig.h"
#include "argcomp.h"


template<class ARFLOAT>
class ARluCompGenEig:
  public virtual
    ARCompGenEig<ARFLOAT, ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
                 ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > > {

 private:

 // a) Data structure used to store matrices.

  ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > Pencil;

 // b) Protected functions:

  virtual void Copy(const ARluCompGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<ARFLOAT> sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<ARFLOAT> sigmap);

 // c.2) Constructors and destructor.

  ARluCompGenEig() { }
  // Short constructor.

  ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
                 const std::string& whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
                 arcomplex<ARFLOAT> sigma, const std::string& whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluCompGenEig(const ARluCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluCompGenEig() { }

 // d) Operators.

  ARluCompGenEig& operator=(const ARluCompGenEig& other);
  // Assignment operator.

}; // class ARluCompGenEig.


// ------------------------------------------------------------------------ //
// ARluCompGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluCompGenEig<ARFLOAT>::
Copy(const ARluCompGenEig<ARFLOAT>& other)
{

  ARCompGenEig<ARFLOAT, ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
               ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;
  if (this->mode > 2) this->objOP->FactorAsB(this->sigmaR);

} // Copy.


template<class ARFLOAT>
inline void ARluCompGenEig<ARFLOAT>::ChangeShift(arcomplex<ARFLOAT> sigmap)
{

  this->objOP->FactorAsB(sigmap);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmap);

} // ChangeShift.


template<class ARFLOAT>
inline void ARluCompGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetRegularMode(&Pencil,
                   &ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluCompGenEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARCompGenEig<ARFLOAT, ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>,
               ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil,
                       &ARluNonSymPencil<arcomplex<ARFLOAT>,ARFLOAT>::MultInvAsBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARluCompGenEig<ARFLOAT>::
ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B, const std::string& whichp,
               int ncvp, ARFLOAT tolp, int maxitp,
               arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv,
                   &Pencil, 
                   &ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluCompGenEig<ARFLOAT>::
ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
               arcomplex<ARFLOAT> sigmap, const std::string& whichp, int ncvp,
               ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvAsBv,
                   &Pencil, 
                   &ARluNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  SetShiftInvertMode(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARluCompGenEig<ARFLOAT>& ARluCompGenEig<ARFLOAT>::
operator=(const ARluCompGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLGCOMP_H
