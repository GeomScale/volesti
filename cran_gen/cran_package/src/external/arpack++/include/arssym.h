/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARSSym.h.
   Arpack++ class ARSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSSYM_H
#define ARSSYM_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arseig.h"
#include "arrssym.h"


template<class ARFLOAT, class ARFOP>
class ARSymStdEig:
  public virtual ARStdEig<ARFLOAT, ARFLOAT, ARFOP>,
  public virtual ARrcSymStdEig<ARFLOAT> {

 public:

 // a) Constructors and destructor.

  ARSymStdEig() { }
  // Short constructor.

  ARSymStdEig(int np, int nevp, ARFOP* objOPp,
              void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
              const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
              int maxitp = 0, ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARSymStdEig(int np, int nevp, ARFOP* objOPp,
              void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
              ARFLOAT sigmap, const std::string& whichp = "LM", int ncvp = 0,
              ARFLOAT tolp = 0.0, int maxitp = 0, ARFLOAT* residp = NULL,
              bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARSymStdEig(const ARSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARSymStdEig& operator=(const ARSymStdEig& other);
  // Assignment operator.

}; // class ARSymStdEig.


// ------------------------------------------------------------------------ //
// ARSymStdEig member functions definition.                                 //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARFOP>
inline ARSymStdEig<ARFLOAT, ARFOP>::
ARSymStdEig(int np, int nevp, ARFOP* objOPp,
            void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
            const std::string& whichp, int ncvp, ARFLOAT tolp,
            int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT, class ARFOP>
inline ARSymStdEig<ARFLOAT, ARFOP>::
ARSymStdEig(int np, int nevp, ARFOP* objOPp,
            void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
            ARFLOAT sigmap, const std::string& whichp, int ncvp, ARFLOAT tolp,
            int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->ChangeShift(sigmap);
  this->DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT, class ARFOP>
ARSymStdEig<ARFLOAT, ARFOP>& ARSymStdEig<ARFLOAT, ARFOP>::
operator=(const ARSymStdEig<ARFLOAT, ARFOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSSYM_H

