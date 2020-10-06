/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARGComp.h.
   Arpack++ class ARCompGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARGCOMP_H
#define ARGCOMP_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arscomp.h"
#include "argeig.h"

template<class ARFLOAT, class ARFOP, class ARFB>
class ARCompGenEig:
  virtual public ARGenEig<ARFLOAT, arcomplex<ARFLOAT>, ARFOP, ARFB>,
  virtual public ARCompStdEig<ARFLOAT, ARFOP>  {

 public:

  // a) Constructors and destructor.

  ARCompGenEig() { }
  // Short constructor (Does nothing but calling base classes constructors).

  ARCompGenEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
               ARFB* objBp,
               void (ARFB::* MultBxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
               const std::string& whichp = "LM", int ncvp = 0,
               ARFLOAT tolp = 0.0, int maxitp = 0,
               arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARCompGenEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
               ARFB* objBp,
               void (ARFB::* MultBxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
               arcomplex<ARFLOAT> sigmap,
               const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
               int maxitp = 0, arcomplex<ARFLOAT>* residp = NULL,
               bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARCompGenEig(const ARCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARCompGenEig() { }
  // Destructor.

 // b) Operators.

  ARCompGenEig& operator=(const ARCompGenEig& other);
  // Assignment operator.

}; // class ARCompGenEig.


// ------------------------------------------------------------------------ //
// ARCompGenEig member functions definition.                                //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARFOP, class ARFB>
inline ARCompGenEig<ARFLOAT, ARFOP, ARFB>::
ARCompGenEig(int np, int nevp, ARFOP* objOPp,
             void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
             ARFB* objBp,
             void (ARFB::* MultBxp)(arcomplex<ARFLOAT>[], arcomplex<ARFLOAT>[]),
             const std::string& whichp, int ncvp, ARFLOAT tolp,
             int maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT, class ARFOP, class ARFB>
inline ARCompGenEig<ARFLOAT, ARFOP, ARFB>::
ARCompGenEig(int np, int nevp, ARFOP* objOPp,
             void (ARFOP::* MultOPxp)(arcomplex<ARFLOAT>[],arcomplex<ARFLOAT>[]),
             ARFB* objBp,
             void (ARFB::* MultBxp)(arcomplex<ARFLOAT>[], arcomplex<ARFLOAT>[]),
             arcomplex<ARFLOAT> sigmap, const std::string& whichp, int ncvp, ARFLOAT tolp,
             int maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->ChangeShift(sigmap);
  this->DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT, class ARFOP, class ARFB>
ARCompGenEig<ARFLOAT, ARFOP, ARFB>& ARCompGenEig<ARFLOAT, ARFOP, ARFB>::
operator=(const ARCompGenEig<ARFLOAT, ARFOP, ARFB>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARGCOMP_H
