/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARSEig.h.
   Arpack++ class ARStdEig definition.
   This class is the base class for all
   standard and generalized problem templates.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSEIG_H
#define ARSEIG_H

#include <cstddef>
#include "arch.h"
#include "arerror.h"
#include "arrseig.h"

// ARStdEig class definition.

template<class ARFLOAT, class ARTYPE, class ARFOP>
class ARStdEig: virtual public ARrcStdEig<ARFLOAT, ARTYPE> {

 public:

 // a) Notation.

  typedef void (ARFOP::* TypeOPx)(ARTYPE[], ARTYPE[]);


 protected:

 // b) User defined parameters.

  ARFOP   *objOP;     // Object that has MultOPx as a member function.
  TypeOPx MultOPx;    // Function that evaluates the product OP*x.

 // c) Protected functions.

  virtual void Copy(const ARStdEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // d) Public functions:

 // d.1) Function that stores user defined parameters.

  virtual void DefineParameters(int np, int nevp, ARFOP* objOPp,
                                TypeOPx MultOPxp, const std::string& whichp="LM",
                                int ncvp=0, ARFLOAT tolp=0.0, int maxitp=0,
                                ARTYPE* residp=NULL, bool ishiftp=true);
  // Set values of problem parameters (also called by constructors).
  // Redefined in ARGenEigenProblem.

 // d.2) Function that allow changes in problem parameters.

  void ChangeMultOPx(ARFOP* objOPp, TypeOPx MultOPxp);
  // Changes the matrix-vector function that performs OP*x.

  virtual void SetRegularMode(ARFOP* objOPp, TypeOPx MultOPxp);
  // Turns problem to regular mode.

  virtual void SetShiftInvertMode(ARTYPE sigmap, ARFOP* objOPp, 
                                  TypeOPx MultOPxp);
  // Turns problem to shift and invert mode with shift defined by sigmap.

 // d.3) Function that permits step by step execution of ARPACK.

  virtual void Iterate() {
    throw ArpackError(ArpackError::NOT_IMPLEMENTED, "Iterate");
  }
  // Takes one iteration of IRA method.


 // d.4) Function that performs all calculations in one step.

  virtual int FindArnoldiBasis();
  // Determines the Arnoldi basis related to the given problem.
  // Redefined in ARGenEigenProblem and ARSymGenEigenProblem.


 // d.5) Constructor and destructor.

  ARStdEig() { }
  // Constructor that does nothing but calling base class constructor.

  ARStdEig(const ARStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARStdEig() { }
  // Very simple destructor.

 // e) Operators.

  ARStdEig& operator=(const ARStdEig& other);
  // Assignment operator.

}; // class ARStdEig.


// ------------------------------------------------------------------------ //
// ARStdEig member functions definition.                                    //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARTYPE, class ARFOP>
inline void ARStdEig<ARFLOAT, ARTYPE, ARFOP>::
Copy(const ARStdEig<ARFLOAT, ARTYPE, ARFOP>& other)
{

  ARrcStdEig<ARFLOAT, ARTYPE>::Copy(other);
  objOP   = other.objOP;
  MultOPx = other.MultOPx;

} // Copy.


template<class ARFLOAT, class ARTYPE, class ARFOP>
void ARStdEig<ARFLOAT, ARTYPE, ARFOP>::
DefineParameters(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARTYPE[], ARTYPE[]), const std::string& whichp,
                 int ncvp, ARFLOAT tolp, int maxitp, ARTYPE* residp, 
                 bool ishiftp)


{

  ARrcStdEig<ARFLOAT, ARTYPE>::DefineParameters(np, nevp, whichp, ncvp, tolp,
                                                maxitp, residp, ishiftp);
  objOP     = objOPp;
  MultOPx   = MultOPxp;

} // DefineParameters.


template<class ARFLOAT, class ARTYPE, class ARFOP>
inline void ARStdEig<ARFLOAT, ARTYPE, ARFOP>::
ChangeMultOPx(ARFOP* objOPp, void (ARFOP::* MultOPxp)(ARTYPE[], ARTYPE[]))
{

  objOP   = objOPp;
  MultOPx = MultOPxp;
  this->Restart();

} // ChangeMultOPx.


template<class ARFLOAT, class ARTYPE, class ARFOP>
inline void ARStdEig<ARFLOAT, ARTYPE, ARFOP>::
SetRegularMode(ARFOP* objOPp, void (ARFOP::* MultOPxp)(ARTYPE[], ARTYPE[]))
{

  ChangeMultOPx(objOPp, MultOPxp);
  this->NoShift();

} // SetRegularMode.


template<class ARFLOAT, class ARTYPE, class ARFOP>
inline void ARStdEig<ARFLOAT, ARTYPE, ARFOP>::
SetShiftInvertMode(ARTYPE sigmap, ARFOP* objOPp, 
                   void (ARFOP::* MultOPxp)(ARTYPE[], ARTYPE[]))
{

  ChangeMultOPx(objOPp, MultOPxp);
  this->ChangeShift(sigmap);

} // SetShiftInvertMode.


template<class ARFLOAT, class ARTYPE, class ARFOP>
int ARStdEig<ARFLOAT, ARTYPE, ARFOP>::FindArnoldiBasis()
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

    if ((this->ido == -1) || (this->ido == 1)) {

      // Performing Matrix vector multiplication: y <- OP*x.

      (objOP->*MultOPx)(&this->workd[this->ipntr[1]],&this->workd[this->ipntr[2]]);

    }

  }
  return this->nconv;

} // FindArnoldiBasis.


template<class ARFLOAT, class ARTYPE, class ARFOP>
ARStdEig<ARFLOAT, ARTYPE, ARFOP>& ARStdEig<ARFLOAT, ARTYPE, ARFOP>::
operator=(const ARStdEig<ARFLOAT, ARTYPE, ARFOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSEIG_H

