/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARRGComp.h.
   Arpack++ class ARrcCompGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRGCOMP_H
#define ARRGCOMP_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arrscomp.h"
#include "arrgeig.h"

template<class ARFLOAT>
class ARrcCompGenEig:
  virtual public ARrcGenEig<ARFLOAT, arcomplex<ARFLOAT> >,
  virtual public ARrcCompStdEig<ARFLOAT>  {

 public:

  // a) Constructors and destructor.

  ARrcCompGenEig() { }
  // Short constructor (Does nothing but calling base classes constructors).

  ARrcCompGenEig(int np, int nevp, const std::string& whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcCompGenEig(int np, int nevp, arcomplex<ARFLOAT> sigmap,
                 const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
                 int maxitp = 0, arcomplex<ARFLOAT>* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcCompGenEig(const ARrcCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcCompGenEig() { }
  // Destructor.

 // b) Operators.

  ARrcCompGenEig& operator=(const ARrcCompGenEig& other);
  // Assignment operator.

}; // class ARrcCompGenEig.


// ------------------------------------------------------------------------ //
// ARrcCompGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline ARrcCompGenEig<ARFLOAT>::
ARrcCompGenEig(int np, int nevp, const std::string& whichp, int ncvp, ARFLOAT tolp,
               int maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARrcCompGenEig<ARFLOAT>::
ARrcCompGenEig(int np, int nevp, arcomplex<ARFLOAT> sigmap, const std::string& whichp,
               int ncvp, ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  this->ChangeShift(sigmap);
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shif and invert mode).


template<class ARFLOAT>
ARrcCompGenEig<ARFLOAT>& ARrcCompGenEig<ARFLOAT>::
operator=(const ARrcCompGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRGCOMP_H

