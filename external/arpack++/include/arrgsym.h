/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARRGSym.h.
   Arpack++ class ARrcSymGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRGSYM_H
#define ARRGSYM_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arrssym.h"
#include "arrgeig.h"

template<class ARFLOAT>
class ARrcSymGenEig:
  virtual public ARrcGenEig<ARFLOAT, ARFLOAT>,
  virtual public ARrcSymStdEig<ARFLOAT> {

 protected:

 // a) Protected variable:

  char    InvertMode;


 // b) Protected functions:

  char CheckInvertMode(char InvertModep);

  virtual void Copy(const ARrcSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  void ChangeInvertMode(char InvertModep);
  // Changes "InvertMode" to 'S' (shift-and-invert),
  // 'B' (buckling) or 'C' (cayley) mode.

  virtual void ChangeShift(ARFLOAT sigmap);
  // Changes shift value.

  virtual void SetShiftInvertMode(ARFLOAT sigmap);
  // Turns problem to shift and invert mode with shift defined by sigmap.

  virtual void SetBucklingMode(ARFLOAT sigmap);
  // Turns problem to buckling mode with shift defined by sigmap.

  virtual void SetCayleyMode(ARFLOAT sigmap);
  // Turns problem to Cayley mode with shift defined by sigmap.


 // c.2) Constructors and destructor.

  ARrcSymGenEig() { InvertMode = 'S'; }
  // Short constructor that does almost nothing.

  ARrcSymGenEig(int np, int nevp, const std::string& whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcSymGenEig(char invertmodep, int np, int nevp, ARFLOAT sigmap,
                const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
                int maxitp = 0, ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift-and-invert, buckling and Cayley modes).

  ARrcSymGenEig(const ARrcSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARrcSymGenEig& operator=(const ARrcSymGenEig& other);
  // Assignment operator.

}; // class ARrcSymGenEig.


// ------------------------------------------------------------------------ //
// ARrcSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline char ARrcSymGenEig<ARFLOAT>::CheckInvertMode(char InvertModep)
{
  if ((InvertModep != 'S') && (InvertModep != 'B') && (InvertModep != 'C')) {
    throw ArpackError(ArpackError::INVMODE_UNDEFINED);
  }
  return InvertModep;

} // CheckInvertMode.


template<class ARFLOAT>
inline void ARrcSymGenEig<ARFLOAT>::Copy(const ARrcSymGenEig<ARFLOAT>& other)
{

  ARrcStdEig<ARFLOAT, ARFLOAT>::Copy(other);
  InvertMode = other.InvertMode;

} // Copy.


template<class ARFLOAT>
inline void ARrcSymGenEig<ARFLOAT>::ChangeInvertMode(char InvertModep)
{

  InvertMode = CheckInvertMode(InvertModep);
  switch (InvertMode) {
  case 'S':
    this->mode    = 3;    // Shift and invert mode.
    break;
  case 'B':
    this->mode    = 4;    // Buckling mode.
    break;
  case 'C':
    this->mode    = 5;    // Cayley mode.
    break;
  }
  this->iparam[7] = this->mode;
  this->Restart();

} // ChangeInvertMode.


template<class ARFLOAT>
inline void ARrcSymGenEig<ARFLOAT>::ChangeShift(ARFLOAT sigmap)
{

  this->sigmaR    = sigmap;
  this->sigmaI    = 0.0;
  ChangeInvertMode(InvertMode);

} // ChangeShift.


template<class ARFLOAT>
void ARrcSymGenEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)

{

  InvertMode = 'S';
  ChangeShift(sigmap);

} // SetShiftInvertMode.


template<class ARFLOAT>
void ARrcSymGenEig<ARFLOAT>::SetBucklingMode(ARFLOAT sigmap)

{

  InvertMode = 'B';
  ChangeShift(sigmap);

} // SetBucklingMode.


template<class ARFLOAT>
void ARrcSymGenEig<ARFLOAT>::SetCayleyMode(ARFLOAT sigmap)

{

  InvertMode = 'C';
  ChangeShift(sigmap);

} // SetCayleyMode.


template<class ARFLOAT>
inline ARrcSymGenEig<ARFLOAT>::
ARrcSymGenEig(int np, int nevp, const std::string& whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  InvertMode = 'S';   // Considering mode = 3 in ChangeShift.
  this->NoShift();
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARrcSymGenEig<ARFLOAT>::
ARrcSymGenEig(char InvertModep, int np, int nevp,
              ARFLOAT sigmap, const std::string& whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  InvertMode = CheckInvertMode(InvertModep); // InvertMode = 'S', 'B', 'C'.
  ChangeShift(sigmap);
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift-and-invert, buckling and Cayley modes).


template<class ARFLOAT>
ARrcSymGenEig<ARFLOAT>& ARrcSymGenEig<ARFLOAT>::
operator=(const ARrcSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRGSYM_H

