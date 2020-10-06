/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARSNSym.h.
   Arpack++ class ARNonSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSNSYM_H
#define ARSNSYM_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arseig.h"
#include "arrsnsym.h"


template<class ARFLOAT, class ARFOP>
class ARNonSymStdEig:
  public virtual ARStdEig<ARFLOAT, ARFLOAT, ARFOP>,
  public virtual ARrcNonSymStdEig<ARFLOAT> {

 public:

 // a) Constructors and destructor.

  ARNonSymStdEig() { }
  // Short constructor.

  ARNonSymStdEig(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
                 const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
                 int maxitp = 0, ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARNonSymStdEig(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
                 ARFLOAT sigma, const std::string& whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0, ARFLOAT* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARNonSymStdEig(const ARNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARNonSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARNonSymStdEig& operator=(const ARNonSymStdEig& other);
  // Assignment operator.

}; // class ARNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARNonSymStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARFOP>
inline ARNonSymStdEig<ARFLOAT, ARFOP>::
ARNonSymStdEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
               const std::string& whichp, int ncvp, ARFLOAT tolp, int maxitp,
               ARFLOAT* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT, class ARFOP>
inline ARNonSymStdEig<ARFLOAT, ARFOP>::
ARNonSymStdEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
               ARFLOAT sigmap, const std::string& whichp, int ncvp, ARFLOAT tolp,
               int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->ChangeShift(sigmap);
  this->DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT, class ARFOP>
ARNonSymStdEig<ARFLOAT, ARFOP>& ARNonSymStdEig<ARFLOAT, ARFOP>::
operator=(const ARNonSymStdEig<ARFLOAT, ARFOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSNSYM_H
