/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARGEig.h.
   Arpack++ class ARGenEig definition.
   Derived from ARStdEig, this class is the
   base class for all generalized eigenvalue problems definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARGEIG_H
#define ARGEIG_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arerror.h"
#include "arrgeig.h"
#include "arseig.h"

// ARGenEig class definition.

template<class ARFLOAT, class ARTYPE, class ARFOP, class ARFB>
class ARGenEig:
  virtual public ARrcGenEig<ARFLOAT, ARTYPE>,
  virtual public ARStdEig<ARFLOAT, ARTYPE, ARFOP> {

 public:

 // a) Notation.

  typedef void (ARFB::* TypeBx)(ARTYPE[], ARTYPE[]);
  typedef void (ARFOP::* TypeOPx)(ARTYPE[], ARTYPE[]);


 protected:

 // b) Protected variables:

  ARFB    *objB;      // Object that has MultBx as a member function.
  TypeBx  MultBx;     // Function that evaluates the product B*x.

 // c) Protected functions:

  virtual void Copy(const ARGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // d) Public functions:

 // d.1) Function that stores user defined parameters.

  virtual void DefineParameters(int np, int nevp, ARFOP* objOPp,
                                TypeOPx MultOPxp, ARFB* objBp, 
                                TypeBx MultBxp, const std::string& whichp="LM", 
                                int ncvp=0, ARFLOAT tolp=0.0,
                                int maxitp=0, ARTYPE* residp=NULL,
                                bool ishiftp=true);
  // Set values of problem parameters (also called by constructors).


 // d.2) Function that allow changes in problem parameters.

  void ChangeMultBx(ARFB* objBp, TypeBx MultBxp);
  // Changes the matrix-vector function that performs B*x.


 // d.3) Functions that perform all calculations in one step.

  virtual int FindArnoldiBasis();
  // Determines the Arnoldi basis related to the given problem.


 // d.4) Constructors and destructor.

  ARGenEig() { }
  // Constructor that does nothing but calling base classes constructors.

  ARGenEig(const ARGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARGenEig() { }
  // Destructor (presently meaningless).

 // e) Operators.

  ARGenEig& operator=(const ARGenEig& other);
  // Assignment operator.

}; // class ARGenEig.


// ------------------------------------------------------------------------ //
// ARGenEig member functions definition.                                    //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARTYPE, class ARFOP, class ARFB>
inline void ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>::
Copy(const ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>& other)
{

  ARStdEig<ARFLOAT, ARTYPE, ARFOP>::Copy(other);
  objB   = other.objB;
  MultBx = other.MultBx;

} // Copy.


template<class ARFLOAT, class ARTYPE, class ARFOP, class ARFB>
void ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>::
DefineParameters(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARTYPE[], ARTYPE[]), ARFB* objBp,
                 void (ARFB::* MultBxp)(ARTYPE[], ARTYPE[]), const std::string& whichp,
                 int ncvp, ARFLOAT tolp, int maxitp, ARTYPE* residp, 
                 bool ishiftp)

{

  // Setting parameters of generalized problems.

  objB   = objBp;
  MultBx = MultBxp;

  // Setting common eigen-problem parameters.

  ARStdEig<ARFLOAT, ARTYPE, ARFOP>::
    DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                     ncvp, tolp, maxitp, residp, ishiftp);

} // DefineParameters.


template<class ARFLOAT, class ARTYPE, class ARFOP, class ARFB>
inline void ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>::
ChangeMultBx(ARFB* objBp, void (ARFB::* MultBxp)(ARTYPE[], ARTYPE[]))
{

  objB   = objBp;
  MultBx = MultBxp;
  this->Restart();

} // ChangeMultBx.


template<class ARFLOAT, class ARTYPE, class ARFOP, class ARFB>
int ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>::FindArnoldiBasis()
{

  if (!this->BasisOK) this->Restart();

  // Changing to auto shift mode.

  if (!this->AutoShift) {
    ArpackError::Set(ArpackError::CHANGING_AUTOSHIFT, "FindArnoldiBasis");
    this->AutoShift=true;
  }

  // ARPACK main loop.

  while (!this->BasisOK) {

    // Calling Aupp.

    try { this->TakeStep(); }
    catch (ArpackError) {
      ArpackError(ArpackError::CANNOT_FIND_BASIS, "FindArnoldiBasis");
      return 0;
    }

    switch (this->ido) {
    case -1:

      // Performing y <- OP*B*x for the first time when mode != 2.

      if (this->mode != 2) {
        this->ipntr[3] = this->ipntr[2]+this->n; // not a clever idea, but...
        (this->objB->*MultBx)(&this->workd[this->ipntr[1]],&this->workd[this->ipntr[3]]);
      }

    case  1:

      // Performing y <- OP*w.

      if (this->mode == 2) { // w = x if mode = 2.
        (this->objOP->*(this->MultOPx))(&this->workd[this->ipntr[1]],&this->workd[this->ipntr[2]]);
      }
      else {           // w = B*x otherwise.
        (this->objOP->*(this->MultOPx))(&this->workd[this->ipntr[3]],&this->workd[this->ipntr[2]]);
      }
      break;

    case  2:

      // Performing y <- B*x.

      (this->objB->*MultBx)(&this->workd[this->ipntr[1]],&this->workd[this->ipntr[2]]);

    }
  }
  return this->nconv;

} // FindArnoldiBasis.


template<class ARFLOAT, class ARTYPE, class ARFOP, class ARFB>
ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>& ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>::
operator=(const ARGenEig<ARFLOAT, ARTYPE, ARFOP, ARFB>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARGEIG_H

