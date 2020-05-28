/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARRGNSym.h.
   Arpack++ class ARrcNonSymGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRGNSYM_H
#define ARRGNSYM_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arrsnsym.h"
#include "arrgeig.h"

template<class ARFLOAT>
class ARrcNonSymGenEig:
  virtual public ARrcGenEig<ARFLOAT, ARFLOAT>,
  virtual public ARrcNonSymStdEig<ARFLOAT>  {

 protected:

 // a) Protected variables:

  char part;


 // b) Protected functions:

  char CheckPart(char partp);

  virtual void Copy(const ARrcNonSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that provides access to internal variables' values.

  ARFLOAT GetShiftImag() { return this->sigmaI; }
  // Returns the imaginary part of the shift (when in shift and invert mode).


 // c.2) Functions that allow changes in problem parameters.

  void ChangePart(char partp);
  // Changes "part" to 'R' (real) or 'I' (imaginary).

  virtual void ChangeShift(ARFLOAT sigmaRp, ARFLOAT sigmaIp = 0.0);
  // Turns the problem to shift-and-invert mode
  // with shift defined by sigmaRp and sigmaIp.

  virtual void SetShiftInvertMode(ARFLOAT sigmaRp);
  // Turns the problem to real shift-and-invert mode with sigmaRp as shift.

  virtual void SetComplexShiftMode(char partp, ARFLOAT sigmaRp, 
                                   ARFLOAT sigmaIp);
  // Turns the problem to complex shift-and-invert mode with shift
  // defined by sigmaRp and sigmaIp.


 // c.3) Constructors and destructor.

  ARrcNonSymGenEig() { part = 'R'; }
  // Short constructor that does almost nothing.

  ARrcNonSymGenEig(int np, int nevp, const std::string& whichp = "LM",
                   int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcNonSymGenEig(int np, int nevp, ARFLOAT sigmap,
                   const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
                   int maxitp = 0, ARFLOAT* residp = NULL, 
                   bool ishiftp = true);
  // Long constructor (real shift and invert mode).

  ARrcNonSymGenEig(int np, int nevp,
                   char partp, ARFLOAT sigmaRp, ARFLOAT sigmaIp,
                   const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
                   int maxitp = 0, ARFLOAT* residp = NULL, 
                   bool ishiftp = true);
  // Long constructor (complex shift and invert mode).

  ARrcNonSymGenEig(const ARrcNonSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcNonSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARrcNonSymGenEig& operator=(const ARrcNonSymGenEig& other);
  // Assignment operator.

}; // class ARrcNonSymGenEig.


// ------------------------------------------------------------------------ //
// ARrcNonSymGenEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline char ARrcNonSymGenEig<ARFLOAT>::CheckPart(char partp)
{
  if ((partp != 'R') && (partp != 'I')) {
    throw ArpackError(ArpackError::PART_UNDEFINED);
  }
  return partp;
} // CheckPart.


template<class ARFLOAT>
inline void ARrcNonSymGenEig<ARFLOAT>::
Copy(const ARrcNonSymGenEig<ARFLOAT>& other)
{

  ARrcStdEig<ARFLOAT, ARFLOAT>::Copy(other);
  part = other.part;

} // Copy.


template<class ARFLOAT>
inline void ARrcNonSymGenEig<ARFLOAT>::ChangePart(char partp)
{

  part = CheckPart(partp);
  if (part == 'R') {
    this->mode    = 3;    // Real part.
  }
  else {
    this->mode    = 4;    // Imaginary part.
  }
  this->iparam[7] = this->mode;
  this->Restart();

} // ChangePart.


template<class ARFLOAT>
inline void ARrcNonSymGenEig<ARFLOAT>::
ChangeShift(ARFLOAT sigmaRp, ARFLOAT sigmaIp)
{

  this->sigmaR    = sigmaRp;
  this->sigmaI    = sigmaIp;
  ChangePart(part);

} // ChangeShift.


template<class ARFLOAT>
inline void ARrcNonSymGenEig<ARFLOAT>::
SetShiftInvertMode(ARFLOAT sigmaRp)
{

  part = 'R';
  ChangeShift(sigmaRp);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline void ARrcNonSymGenEig<ARFLOAT>::
SetComplexShiftMode(char partp, ARFLOAT sigmaRp, ARFLOAT sigmaIp)
{

  part   = CheckPart(partp);
  ChangeShift(sigmaRp, sigmaIp);

} // SetComplexShiftMode.


template<class ARFLOAT>
inline ARrcNonSymGenEig<ARFLOAT>::
ARrcNonSymGenEig(int np, int nevp, const std::string& whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)
{

  part = 'R';                // Considering mode = 3 in ChangeShift.
  this->NoShift();
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARrcNonSymGenEig<ARFLOAT>::
ARrcNonSymGenEig(int np, int nevp, ARFLOAT sigmap, const std::string& whichp, int ncvp,
                 ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)
{

  SetShiftInvertMode(sigmap);
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);


} // Long constructor (real shift and invert mode).


template<class ARFLOAT>
inline ARrcNonSymGenEig<ARFLOAT>::
ARrcNonSymGenEig(int np, int nevp, char partp, ARFLOAT sigmaRp,
                 ARFLOAT sigmaIp, const std::string& whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)
{

  SetComplexShiftMode(partp, sigmaRp, sigmaIp);
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARrcNonSymGenEig<ARFLOAT>& ARrcNonSymGenEig<ARFLOAT>::
operator=(const ARrcNonSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRGNSYM_H

