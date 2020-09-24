/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARSComp.h.
   Arpack++ class ARCompStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSCOMP_H
#define ARSCOMP_H

#include <cstddef>

#include "arch.h"
#include "arseig.h"
#include "arrscomp.h"

template<class ARFLOAT, class ARFOP>
class ARCompStdEig:
  virtual public ARStdEig<ARFLOAT, arcomplex<ARFLOAT>, ARFOP>,
  virtual public ARrcCompStdEig<ARFLOAT> {

 public:

 // a) Constructors and destructor.

  ARCompStdEig() { }
  // Short constructor.

  ARCompStdEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
               const std::string& whichp = "LM", int ncvp = 0,
               ARFLOAT tolp = 0.0, int maxitp = 0,
               arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARCompStdEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
               arcomplex<ARFLOAT> sigma,  const std::string& whichp = "LM",
               int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
               arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARCompStdEig(const ARCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARCompStdEig() { }
  // Destructor.

 // b) Operators.

  ARCompStdEig& operator=(const ARCompStdEig& other);
  // Assignment operator.

}; // class ARCompStdEig.


// ------------------------------------------------------------------------ //
// ARCompStdEig member functions definition.                                //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARFOP>
inline ARCompStdEig<ARFLOAT, ARFOP>::
ARCompStdEig(int np, int nevp, ARFOP* objOPp,
             void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
             const std::string& whichp, int ncvp, ARFLOAT tolp, int maxitp,
             arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT, class ARFOP>
inline ARCompStdEig<ARFLOAT, ARFOP>::
ARCompStdEig(int np, int nevp, ARFOP* objOPp,
             void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
             arcomplex<ARFLOAT> sigmap, const std::string& whichp, int ncvp,
             ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
             bool ishiftp)

{

  this->ChangeShift(sigmap);
  this->DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT, class ARFOP>
ARCompStdEig<ARFLOAT, ARFOP>& ARCompStdEig<ARFLOAT, ARFOP>::
operator=(const ARCompStdEig<ARFLOAT, ARFOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSCOMP_H
