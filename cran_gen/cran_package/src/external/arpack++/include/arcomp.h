/*
  ARPACK++ v1.2 2/20/2000
  c++ interface to ARPACK code.

  MODULE arcomp.h
  arcomplex complex type definition.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef ARCOMP_H
#define ARCOMP_H

#include <complex>

#ifdef __GNUG__
  
#define arcomplex std::complex

#endif

#if defined(__SUNPRO_CC) || defined(__sgi)

  template <class ARFLOAT>
  class arcomplex: public complex
  {
   public:

    arcomplex(ARFLOAT x, ARFLOAT y): complex(x,y) { }
    arcomplex(): complex() { }
    arcomplex(complex x): complex(x) { }

  };

#endif

#endif // ARCOMP_H



