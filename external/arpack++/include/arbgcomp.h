/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBGComp.h.
   Arpack++ class ARluCompGenEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBGCOMP_H
#define ARBGCOMP_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arbnsmat.h"
#include "arbnspen.h"
#include "arrseig.h"
#include "argcomp.h"


template<class ARFLOAT>
class ARluCompGenEig:
  public virtual
    ARCompGenEig<ARFLOAT, ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
                 ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > > {

 private:

 // a) Data structure used to store matrices.

  ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > Pencil;

 // b) Protected functions:

  virtual void Copy(const ARluCompGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<ARFLOAT> sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<ARFLOAT> sigmap);

 // c.2) Constructors and destructor.

  ARluCompGenEig() { }
  // Short constructor.

  ARluCompGenEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
                 const std::string& whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluCompGenEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
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

  ARCompGenEig<ARFLOAT, ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
               ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARluCompGenEig<ARFLOAT>::
ChangeShift(arcomplex<ARFLOAT> sigmaRp)
{

  this->objOP->FactorAsB(sigmaRp);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmaRp);

} // ChangeShift.


template<class ARFLOAT>
inline void ARluCompGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetRegularMode(&Pencil,
                   &ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluCompGenEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARCompGenEig<ARFLOAT, ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>,
               ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil,
                       &ARbdNonSymPencil<arcomplex<ARFLOAT>,ARFLOAT>::MultInvAsBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARluCompGenEig<ARFLOAT>::
ARluCompGenEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B, const std::string& whichp,
               int ncvp, ARFLOAT tolp, int maxitp,
               arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv,
                   &Pencil, 
                   &ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluCompGenEig<ARFLOAT>::
ARluCompGenEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
               arcomplex<ARFLOAT> sigmap, const std::string& whichp, int ncvp,
               ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvAsBv,
                   &Pencil, 
                   &ARbdNonSymPencil<arcomplex<ARFLOAT>,ARFLOAT>::MultBv, 
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


#endif // ARBGCOMP_H
