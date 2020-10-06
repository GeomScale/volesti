/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARRGEig.h.
   Arpack++ class ARrcGenEig definition.
   Derived from ARrcStdEig, this class implements the
   reverse communication interface for generalized problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRGEIG_H
#define ARRGEIG_H

#include "arch.h"
#include "arerror.h"
#include "arrseig.h"

// ARrcGenEig class definition.

template<class ARFLOAT, class ARTYPE>
class ARrcGenEig: virtual public ARrcStdEig<ARFLOAT, ARTYPE> {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  void NoShift();
  // Turns the problem to regular mode.


 // a.2) Constructors and destructor.

  ARrcGenEig();
  // Short constructor that does almost nothing.

  ARrcGenEig(const ARrcGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcGenEig() { }
  // Destructor (presently meaningless).

 // b) Operators.

  ARrcGenEig& operator=(const ARrcGenEig& other);
  // Assignment operator.

}; // class ARrcGenEig.


// ------------------------------------------------------------------------ //
// ARrcGenEig member functions definition.                                  //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARTYPE>
inline void ARrcGenEig<ARFLOAT, ARTYPE>::NoShift()
{

  this->sigmaR    = (ARTYPE)0;
  this->sigmaI    = 0.0;
  this->mode      = 2;
  this->iparam[7] = this->mode;
  this->Restart();

} // NoShift.


template<class ARFLOAT, class ARTYPE>
inline ARrcGenEig<ARFLOAT, ARTYPE>::ARrcGenEig()
{

  this->bmat = 'G';   // This is a generalized problem.
  this->NoShift();

} // Short constructor.


template<class ARFLOAT, class ARTYPE>
ARrcGenEig<ARFLOAT, ARTYPE>& ARrcGenEig<ARFLOAT, ARTYPE>::
operator=(const ARrcGenEig<ARFLOAT, ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRGEIG_H

