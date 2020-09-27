/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLGSym.h.
   Arpack++ class ARluSymGenEig definition
   (SuperLU version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Kristi Maschhoff
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLGSYM_H
#define ARLGSYM_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arlsmat.h"
#include "arlspen.h"
#include "argsym.h"


template<class ARFLOAT>
class ARluSymGenEig:
  public virtual ARSymGenEig<ARFLOAT, ARluSymPencil<ARFLOAT>,
                             ARluSymPencil<ARFLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARluSymPencil<ARFLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARluSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

  virtual void SetBucklingMode(ARFLOAT sigmap);

  virtual void SetCayleyMode(ARFLOAT sigmap);

 // c.2) Constructors and destructor.

  ARluSymGenEig() { }
  // Short constructor.

  ARluSymGenEig(int nevp, ARluSymMatrix<ARFLOAT>& A,
                ARluSymMatrix<ARFLOAT>& B, const std::string& whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluSymGenEig(char InvertModep, int nevp, ARluSymMatrix<ARFLOAT>& A,
                ARluSymMatrix<ARFLOAT>& B, ARFLOAT sigma, const std::string& whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert, buckling and Cayley modes).

  ARluSymGenEig(const ARluSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARluSymGenEig& operator=(const ARluSymGenEig& other);
  // Assignment operator.

}; // class ARluSymGenEig.


// ------------------------------------------------------------------------ //
// ARluSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
Copy(const ARluSymGenEig<ARFLOAT>& other)
{

  ARSymGenEig<ARFLOAT, ARluSymPencil<ARFLOAT>,
              ARluSymPencil<ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;
  this->objA   = &Pencil;
  if (this->mode > 2) this->objOP->FactorAsB(this->sigmaR);

} // Copy.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::ChangeShift(ARFLOAT sigmap)
{

  this->objOP->FactorAsB(sigmap);
  ARrcSymGenEig<ARFLOAT>::ChangeShift(sigmap);

} // ChangeShift.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARluSymPencil<ARFLOAT> >::
    SetRegularMode(&Pencil, &ARluSymPencil<ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
SetShiftInvertMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARluSymPencil<ARFLOAT>, ARluSymPencil<ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARluSymPencil<ARFLOAT>::MultInvAsBv);
  this->ChangeMultBx(&Pencil, &ARluSymPencil<ARFLOAT>::MultBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
SetBucklingMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARluSymPencil<ARFLOAT>, ARluSymPencil<ARFLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARluSymPencil<ARFLOAT>::MultInvAsBv);
  this->ChangeMultBx(&Pencil, &ARluSymPencil<ARFLOAT>::MultAv);

} // SetBucklingMode.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
SetCayleyMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARluSymPencil<ARFLOAT>, ARluSymPencil<ARFLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARluSymPencil<ARFLOAT>::MultInvAsBv,
                  &Pencil, &ARluSymPencil<ARFLOAT>::MultAv);
  this->ChangeMultBx(&Pencil, &ARluSymPencil<ARFLOAT>::MultBv);

} // SetCayleyMode.


template<class ARFLOAT>
inline ARluSymGenEig<ARFLOAT>::
ARluSymGenEig(int nevp, ARluSymMatrix<ARFLOAT>& A,
              ARluSymMatrix<ARFLOAT>& B, const std::string& whichp, int ncvp,
              ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->InvertMode = 'S';
  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARluSymPencil<ARFLOAT>::MultInvBAv, &Pencil,
                   &ARluSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluSymGenEig<ARFLOAT>::
ARluSymGenEig(char InvertModep, int nevp, ARluSymMatrix<ARFLOAT>& A,
              ARluSymMatrix<ARFLOAT>& B, ARFLOAT sigmap,
              const std::string& whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARluSymPencil<ARFLOAT>::MultInvAsBv, &Pencil,
                   &ARluSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  this->InvertMode = this->CheckInvertMode(InvertModep);
  switch (this->InvertMode) {
  case 'B':  // Buckling mode.
    this->ChangeMultBx(&Pencil, &ARluSymPencil<ARFLOAT>::MultAv);
  case 'S':  // Shift and invert mode.
    ChangeShift(sigmap);
    break;
  case 'C':  // Cayley mode.
    SetCayleyMode(sigmap);
  }

} // Long constructor (shift and invert, buckling and Cayley modes).


template<class ARFLOAT>
ARluSymGenEig<ARFLOAT>& ARluSymGenEig<ARFLOAT>::
operator=(const ARluSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLGSYM_H
